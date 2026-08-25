package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Unified long-read adapter trimmer: terminal clipping and chimera splitting
 * in a single step using banded edit distance.
 *
 * <p>Supersedes the separate LONGREADCLIP + LONGREADSPLIT steps from earlier
 * 0.42 iterations.  The key improvements over those steps are:
 * <ol>
 *   <li><b>Edit distance instead of Hamming.</b>  ONT R9 and PacBio CLR errors
 *       are predominantly indels.  A single insertion shifts every subsequent
 *       base out of frame, causing Hamming distance to over-count mismatches
 *       and reject real adapter hits.  Edit distance handles indels correctly.</li>
 *   <li><b>k-mer seed size 6 (was 8).</b>  At 10 % error rate, the probability
 *       that an 8-mer survives error-free is ≈ 0.43; for a 6-mer it is ≈ 0.53.
 *       Shorter seeds improve recall at the cost of more candidates to verify,
 *       which remains fast because the verification step is tightly bounded.</li>
 *   <li><b>Re-clipping after each split.</b>  In the old two-step pipeline,
 *       LONGREADCLIP ran as a separate step and never re-examined the new
 *       termini produced by LONGREADSPLIT.  This left partial adapter sequence
 *       at split junctions.  Here, each split fragment's ends are re-scanned
 *       immediately, removing residual adapter before the fragment is emitted.</li>
 *   <li><b>Platform hint is behaviour-neutral.</b>  {@code ONT}, {@code CLR} and
 *       {@code HIFI} share one path and produce identical output.  The parameter is
 *       parsed and validated, so an unrecognised value is rejected, but nothing
 *       branches on it.  Gating interior splitting by platform was measured against
 *       exact adapter coordinates and costs interior recall entirely, plus 3′ recall
 *       with it, since the near-terminal logic shares the interior scan's code path.
 *       Specificity comes from the seed and full-adapter requirements below, which
 *       hold mid-read over-clipping to 0 and 10 bases across an error range from
 *       0.5% to 15%.</li>
 * </ol>
 *
 * <p><b>Interior vs terminal matching.</b>  Terminal clipping accepts partial
 * adapter overlaps down to {@code minOverlap} bases (the adapter may hang off
 * the read end).  Interior chimera splitting requires the full adapter length:
 * {@code bestMatchAt()} only accepts an overlap equal to {@code adapterLen}.
 * This eliminates false-positive splits caused by short adapter-like motifs in
 * high-error-rate reads, where partial matches (10 bp at 15% = 1 edit) are
 * easily satisfied by chance.
 *
 * <p><b>Zero-allocation k-mer scan.</b>  The interior scan loop converts each
 * read position to a 12-bit integer key (2 bits per base, A=0 C=1 G=2 T=3)
 * and looks up a flat {@code List<int[]>[4096]} table, with no String creation and no
 * HashMap boxing.  For a 50 kb read this eliminates ~50 000 per-position String
 * allocations per thread, substantially reducing GC heap pressure under 40-thread
 * parallelism.
 *
 * <p><b>Single-end mode only.</b>  No current long-read platform produces
 * paired-end data; invoking this step in PE mode throws immediately.
 *
 * <p><b>Known limitation: the edit budget is a hard floor on recall.</b>  A match
 * is accepted only within {@code allowedEdits(len)} edits, which for an 18 bp
 * adapter at {@code maxErrorRate} 0.15 is 2.  An adapter carrying two mismatches
 * of its own plus a sequencing error in the neighbouring bases exceeds that and
 * is not found.  Nothing in this class can recover such a match; it is a property
 * of the threshold, not of where the scans look.
 *
 * <p>Raising the budget was measured as a way to reach those reads and rejected.
 * Rounding instead of truncating (18 bp: 2 edits to 3) cleared almost all of
 * them, but multiplied genuine over-clipping of payload roughly eightfold and
 * lowered overall F1, because a looser threshold admits adapter-like payload just
 * as readily as degraded adapter.  Any unseeded brute-force scan over a window of
 * start positions fails the same way, which is why neither terminal scan uses one.
 *
 * <p>Closing that gap needs a discriminator other than the edit threshold:
 * base qualities at the match, a positional prior, or complexity weighting so
 * that low-complexity matches are required to show more evidence than
 * high-complexity ones.  A uniform threshold cannot separate the two.
 *
 * <p>Recommended pipeline order:
 * <pre>
 *   LONGREADTRIM:adapters.fa:0.10  MINLEN:200
 * </pre>
 *
 * <p>Step syntax:
 * <pre>
 *   LONGREADTRIM:&lt;fasta&gt;:&lt;maxErrorRate&gt;[:&lt;minOverlap&gt;[:&lt;minFragLen&gt;[:&lt;platform&gt;]]]
 * </pre>
 *
 * <p>Parameters:
 * <ul>
 *   <li>{@code fasta}        – FASTA file of adapter sequences. Both orientations
 *       are loaded automatically.</li>
 *   <li>{@code maxErrorRate} – maximum edit-distance fraction allowed in a confirmed
 *       adapter match (e.g. {@code 0.10} = 1 edit per 10 bp).</li>
 *   <li>{@code minOverlap}   – minimum bases of overlap required to call a terminal
 *       match [default 10].  Also defines the terminal exclusion zone for interior
 *       scanning.  Interior matches require the full adapter length.</li>
 *   <li>{@code minFragLen}   – split fragments shorter than this are discarded
 *       [default 100].</li>
 *   <li>{@code platform}     – {@code ONT} (default), {@code CLR}, or {@code HIFI}.
 *       All three behave identically; the value is validated but selects no code
 *       path.  Retained so existing command lines keep working, and because a
 *       future error-rate-dependent behaviour would want it already plumbed.</li>
 * </ul>
 */
public class LongReadTrimmer implements Trimmer {

    // -----------------------------------------------------------------------
    // Constants

    private static final int KMER_SIZE           = 6;
    private static final int KMER_TABLE_SIZE     = 1 << (2 * KMER_SIZE); // 4096
    private static final int DEFAULT_MIN_OVERLAP  = 10;
    private static final int DEFAULT_MIN_FRAG_LEN = 100;

    // -----------------------------------------------------------------------
    // Inner types

    public enum Platform { ONT, CLR, HIFI }

    // -----------------------------------------------------------------------
    // Fields

    /** All adapter orientations (forward + RC). Used for 3′ clipping and interior splitting. */
    private final List<String> adapters;
    /** Forward-orientation adapters only. Used for 5′ clipping. */
    private final List<String> fwdAdapters;
    /**
     * Flat 4096-entry lookup table for 6-mer seeds.  Index is the 12-bit packed
     * encoding of the k-mer (A=0, C=1, G=2, T=3); value is a list of
     * {adapterIndex, positionInAdapter} pairs.  Positions containing N are skipped.
     * Direct array access with integer keys avoids HashMap lookup and autoboxing.
     */
    @SuppressWarnings("unchecked")
    private final List<int[]>[] kmerTable = new List[KMER_TABLE_SIZE];

    private final float    maxErrorRate;
    private final int      minOverlap;
    private final int      minFragLen;
    private final Platform platform;

    /** Per-adapter k-mer presence tables for terminal-scan pre-filtering. */
    private boolean[][] adapterHasKmer;
    private boolean[][] fwdAdapterHasKmer;
    /** True for adapters that contain N wildcards; filter is unsafe for those. */
    private boolean[] adapterHasN;
    private boolean[] fwdAdapterHasN;
    private int maxAdapterLen;
    /** Reusable DP scratch arrays; eliminates per-call int[] allocation for len > 64. */
    private ThreadLocal<int[][]> dpScratch;
    /** Reusable char buffer for the k-mer scans; see readChars(). */
    private ThreadLocal<char[]> charScratch;
    /** Reusable candidate-position buffer for findInternalHits(). */
    private ThreadLocal<int[]> candScratch;

    // -----------------------------------------------------------------------
    // Construction

    public LongReadTrimmer(String args) throws IOException {
        String[] parts = args.split(":");
        if (parts.length < 2)
            throw new IllegalArgumentException(
                    "LONGREADTRIM requires at least <fasta>:<maxErrorRate>");

        Platform tempPlatform = Platform.ONT;
        int      tempMinFrag  = DEFAULT_MIN_FRAG_LEN;
        int      tempMinOvlp  = DEFAULT_MIN_OVERLAP;
        float    tempErrRate;
        int      fileEndIdx;

        int scanIdx = parts.length - 1;

        String rightmost = parts[scanIdx].trim().toUpperCase();
        if (rightmost.equals("ONT") || rightmost.equals("CLR") || rightmost.equals("HIFI")) {
            tempPlatform = Platform.valueOf(rightmost);
            scanIdx--;
        }

        List<Integer> rightInts = new ArrayList<>();
        while (scanIdx >= 0) {
            try {
                rightInts.add(Integer.parseInt(parts[scanIdx].trim()));
                scanIdx--;
            } catch (NumberFormatException e) {
                break;
            }
        }
        switch (rightInts.size()) {
            case 0:  break;
            case 1:  tempMinOvlp = rightInts.get(0); break;
            default: tempMinFrag = rightInts.get(0);
                     tempMinOvlp = rightInts.get(1);
                     break;
        }

        if (scanIdx < 0)
            throw new IllegalArgumentException(
                    "LONGREADTRIM: could not locate maxErrorRate in: " + args);
        try {
            tempErrRate = Float.parseFloat(parts[scanIdx].trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                    "LONGREADTRIM: invalid maxErrorRate '" + parts[scanIdx] + "'");
        }
        fileEndIdx = scanIdx;

        StringBuilder pathBuilder = new StringBuilder();
        for (int i = 0; i < fileEndIdx; i++) {
            if (i > 0) pathBuilder.append(':');
            pathBuilder.append(parts[i]);
        }
        if (pathBuilder.length() == 0)
            throw new IllegalArgumentException("LONGREADTRIM: missing FASTA file path");

        maxErrorRate = tempErrRate;
        minOverlap   = tempMinOvlp;
        minFragLen   = tempMinFrag;
        platform     = tempPlatform;

        adapters    = new ArrayList<>();
        fwdAdapters = new ArrayList<>();
        loadAdapters(new File(pathBuilder.toString()));

        if (adapters.isEmpty())
            throw new IllegalArgumentException(
                    "LONGREADTRIM: no adapter sequences found in " + pathBuilder);

        // Build per-adapter k-mer presence tables for terminal pre-filtering.
        maxAdapterLen = 0;
        for (String a : adapters) maxAdapterLen = Math.max(maxAdapterLen, a.length());

        adapterHasKmer = new boolean[adapters.size()][KMER_TABLE_SIZE];
        for (int ai = 0; ai < adapters.size(); ai++) {
            char[] chars = adapters.get(ai).toCharArray();
            for (int pos = 0; pos <= chars.length - KMER_SIZE; pos++) {
                int code = encodeKmer(chars, pos, KMER_SIZE);
                if (code >= 0) adapterHasKmer[ai][code] = true;
            }
        }
        fwdAdapterHasKmer = new boolean[fwdAdapters.size()][KMER_TABLE_SIZE];
        for (int fi = 0; fi < fwdAdapters.size(); fi++) {
            char[] chars = fwdAdapters.get(fi).toCharArray();
            for (int pos = 0; pos <= chars.length - KMER_SIZE; pos++) {
                int code = encodeKmer(chars, pos, KMER_SIZE);
                if (code >= 0) fwdAdapterHasKmer[fi][code] = true;
            }
        }
        adapterHasN = new boolean[adapters.size()];
        for (int ai = 0; ai < adapters.size(); ai++)
            for (char c : adapters.get(ai).toCharArray())
                if (c == 'N') { adapterHasN[ai] = true; break; }

        fwdAdapterHasN = new boolean[fwdAdapters.size()];
        for (int fi = 0; fi < fwdAdapters.size(); fi++)
            for (char c : fwdAdapters.get(fi).toCharArray())
                if (c == 'N') { fwdAdapterHasN[fi] = true; break; }

        final int scratchLen = maxAdapterLen + 1;
        dpScratch = ThreadLocal.withInitial(() -> new int[][]{ new int[scratchLen], new int[scratchLen] });
        charScratch = ThreadLocal.withInitial(() -> new char[8192]);
        candScratch = ThreadLocal.withInitial(() -> new int[256]);
    }

    // -----------------------------------------------------------------------
    // Adapter loading and k-mer indexing

    private void loadAdapters(File f) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            StringBuilder seq = new StringBuilder();
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.startsWith(">")) {
                    if (seq.length() > 0) { addAdapter(seq.toString()); seq.setLength(0); }
                } else if (!line.isEmpty()) {
                    seq.append(line);
                }
            }
            if (seq.length() > 0) addAdapter(seq.toString());
        }
    }

    private void addAdapter(String seq) {
        String upper = seq.toUpperCase();
        fwdAdapters.add(upper);
        indexAndRegister(upper);
        String rc = reverseComplement(upper);
        if (!rc.equals(upper)) indexAndRegister(rc);
    }

    private void indexAndRegister(String seq) {
        char[] chars = seq.toCharArray();
        int ai = adapters.size();
        adapters.add(seq);
        for (int pos = 0; pos <= chars.length - KMER_SIZE; pos++) {
            int code = encodeKmer(chars, pos, KMER_SIZE);
            if (code < 0) continue; // skip positions with N or ambiguous bases
            if (kmerTable[code] == null) kmerTable[code] = new ArrayList<>();
            kmerTable[code].add(new int[]{ai, pos});
        }
    }

    // -----------------------------------------------------------------------
    // Sequence utilities

    // Package-private so OrientTrimmer can reuse these primitives. Keeping them
    // private would mean a second copy of the same sequence and edit-distance code.
    static String reverseComplement(String seq) {
        int len = seq.length();
        StringBuilder sb = new StringBuilder(len);
        for (int i = len - 1; i >= 0; i--) sb.append(complement(seq.charAt(i)));
        return sb.toString();
    }

    private static char complement(char c) {
        return switch (c) {
            case 'A' -> 'T';
            case 'T' -> 'A';
            case 'G' -> 'C';
            case 'C' -> 'G';
            default  -> 'N';
        };
    }

    /**
     * Encodes {@code k} bases at {@code offset} in {@code chars} as a packed
     * 12-bit integer (A=0, C=1, G=2, T=3).  Returns -1 if any base is
     * ambiguous.  Avoids all heap allocation; the result is used directly as
     * an index into {@link #kmerTable}.
     */
    static int encodeKmer(char[] chars, int offset, int k) {
        int code = 0;
        for (int i = 0; i < k; i++) {
            int b = baseCode(chars[offset + i]);
            if (b < 0) return -1;
            code = (code << 2) | b;
        }
        return code;
    }

    static int baseCode(char c) {
        return switch (c) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default  -> -1;
        };
    }

    // -----------------------------------------------------------------------
    // Edit distance
    //
    // Dispatch: Myers' bit-parallel O(n) for len ≤ 64; DP fallback for len > 64.
    // N bases are wildcards (zero substitution cost in DP; treated as match-all in BP).
    // Returns min(distance, maxEdits+1).

    private int editDistance(String s, int sOff,
                             String t, int tOff,
                             int len, int maxEdits) {
        if (len == 0) return 0;
        if (len <= 64) return editDistanceBP(s, sOff, t, tOff, len, maxEdits);
        return editDistanceDP(s, sOff, t, tOff, len, maxEdits, dpScratch.get());
    }

    /**
     * Myers' bit-parallel edit distance for len ≤ 64.
     * Encodes the pattern (s) into four 64-bit bitmasks (one per base); each text
     * character (t) selects its mask.  N in either string is treated as a wildcard
     * by setting all four bits for that position.
     */
    static int editDistanceBP(String s, int sOff,
                                      String t, int tOff,
                                      int len, int maxEdits) {
        long peqA = 0L, peqC = 0L, peqG = 0L, peqT = 0L;
        for (int j = 0; j < len; j++) {
            long bit = 1L << j;
            switch (s.charAt(sOff + j)) {
                case 'A' -> peqA |= bit;
                case 'C' -> peqC |= bit;
                case 'G' -> peqG |= bit;
                case 'T' -> peqT |= bit;
                default  -> { peqA |= bit; peqC |= bit; peqG |= bit; peqT |= bit; }
            }
        }
        long allOnes = (len == 64) ? -1L : (1L << len) - 1L;
        long Pv = allOnes, Mv = 0L;
        int score = len;
        for (int i = 0; i < len; i++) {
            long Eq = switch (t.charAt(tOff + i)) {
                case 'A' -> peqA;
                case 'C' -> peqC;
                case 'G' -> peqG;
                case 'T' -> peqT;
                default  -> allOnes; // N: wildcard
            };
            long Xv  = Eq | Mv;
            long Ph, Mh;
            long sum  = ((Eq & Pv) + Pv) & allOnes;
            long Xh   = (sum ^ Pv) | Eq;
            Ph   = Mv | (~(Xh | Pv) & allOnes);
            Mh   = Pv & Xh;
            long highBit = 1L << (len - 1);
            if ((Ph & highBit) != 0L) score++;
            if ((Mh & highBit) != 0L) score--;
            Ph = ((Ph << 1) | 1L) & allOnes;
            Mv = Ph & Xv;
            Pv = ((Mh << 1) | ~(Xv | Ph)) & allOnes;
            if (score - (len - 1 - i) > maxEdits) return maxEdits + 1;
        }
        return Math.min(score, maxEdits + 1);
    }

    /** Wagner-Fischer DP fallback for len > 64, reusing caller-supplied scratch arrays. */
    private static int editDistanceDP(String s, int sOff,
                                      String t, int tOff,
                                      int len, int maxEdits,
                                      int[][] scratch) {
        int[] prev = scratch[0];
        int[] curr = scratch[1];
        for (int j = 0; j <= len; j++) prev[j] = j;

        for (int i = 1; i <= len; i++) {
            curr[0] = i;
            int rowMin = i;
            for (int j = 1; j <= len; j++) {
                char sc  = s.charAt(sOff + i - 1);
                char tc  = t.charAt(tOff + j - 1);
                int cost = (sc == tc || sc == 'N' || tc == 'N') ? 0 : 1;
                int val  = Math.min(
                        Math.min(prev[j] + 1,
                                 curr[j - 1] + 1),
                        prev[j - 1] + cost);
                curr[j] = val;
                if (val < rowMin) rowMin = val;
            }
            if (rowMin > maxEdits) return maxEdits + 1;
            int[] tmp = prev; prev = curr; curr = tmp;
        }
        return Math.min(prev[len], maxEdits + 1);
    }

    /** Returns true if any 6-mer from readChars[start..end) exists in adapterKmers. */
    private static boolean hasSharedKmer(char[] readChars, int start, int end,
                                         boolean[] adapterKmers) {
        for (int rp = start; rp <= end - KMER_SIZE; rp++) {
            int code = encodeKmer(readChars, rp, KMER_SIZE);
            if (code < 0) return true; // N in read: wildcard, cannot rule out a match
            if (adapterKmers[code]) return true;
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // Trimmer interface

    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in.length > 1)
            throw new RuntimeException(
                    "LONGREADTRIM cannot be used in paired-end mode. " +
                    "Long-read sequencing platforms (ONT, PacBio) produce single-end data only.");

        if (in[0] == null) return in;

        FastqRecord rec    = in[0];
        String      seq    = rec.getSequence();
        int         seqLen = rec.getLength();

        // ---- Terminal clipping ----
        int clipFrom = findFivePrimeClip(seq, seqLen);
        int trimTo   = findThreePrimeClip(seq, seqLen);

        // ---- Interior adapter hits ----
        // Unconditional: no platform is excluded.  Excluding one costs interior
        // recall entirely, and 3' recall with it, because the near-terminal work in
        // bestMatchAt and findInternalHits is reached only through this path.
        // Measured against exact adapter coordinates at ~99.5% identity, an
        // exclusion here missed 937 of 937 interior adapters and dropped 3' recall
        // to 0.6812.
        //
        // Specificity comes from the seed and full-adapter requirements instead, and
        // it is what makes an unconditional scan safe.  Against LONGREADSPLIT, whose
        // matcher is looser and needs no seed, this path produces 0 bases of mid-read
        // over-clipping where LONGREADSPLIT produces 1352 at ~99.5% identity, and 10
        // against 51014 at 85% identity, both out of ~100 Mbp of payload, at equal
        // recall 1.0000.
        //
        // Scanned on the UNTRIMMED read, in full-read coordinates.  Scanning the
        // already-clipped window made the result depend on the terminal scans
        // running first: a terminal scan can accept a shifted partial alignment
        // and clip through an adapter, and the interior scan then no longer sees
        // a full adapter to remove.  That happens whenever one adapter is nested
        // in another, for LSK114, seq[1..17] of the 18 bp 3' adapter is exactly
        // the first 17 bases of the 5' adapter's reverse complement, so the 3'
        // scan matches the longer adapter one base late and strands the first
        // base.  Scanning the whole read keeps the full-length hit available.
        List<int[]> hits = findInternalHits(seq, seqLen);
        {
            // A hit straddling the 3' boundary wins over the terminal clip: pull
            // trimTo back across the whole adapter instead of keeping the residue.
            //
            // The 3' side only.  The 5' side already has a near-terminal
            // full-adapter scan, so its boundary is trustworthy, and reconciling
            // it there over-clips for the same nesting reason in mirror image:
            // seq[10..26] of the 27 bp 5' adapter equals the first 17 bases of
            // the 3' adapter's reverse complement, so a shifted hit ends one base
            // past the real 5' adapter and would drag clipFrom with it.
            for (int[] hit : hits)
                if (hit[1] >= trimTo && hit[0] < trimTo) trimTo = hit[0];
        }

        int workLen = trimTo - clipFrom;
        if (workLen <= 0) return new FastqRecord[]{null};

        // ---- Interior chimera splitting ----
        // The gate only needs room for the interior scan zone (2 * minOverlap)
        // plus ONE viable fragment: the per-fragment minFragLen filter below
        // discards the non-viable side on its own.  Requiring 2 * minFragLen
        // here skipped the scan on reads that could still yield one clean
        // fragment, leaving the interior adapter in the output.
        if (workLen >= 2 * minOverlap + minFragLen && !hits.isEmpty()) {

            // Only hits fully inside the (possibly widened) kept window split it;
            // the straddling ones were already absorbed into the boundaries.
            List<int[]> inner = new ArrayList<>();
            for (int[] hit : hits)
                if (hit[0] >= clipFrom && hit[1] <= trimTo) inner.add(hit);

            if (!inner.isEmpty()) {
                List<FastqRecord> fragments = new ArrayList<>();
                int start   = 0;                       // work coordinates
                int total   = inner.size() + 1;
                int fragNum = 1;

                for (int[] hit : inner) {
                    int hitStart = hit[0] - clipFrom;   // to work coordinates
                    int fragLen  = hitStart - start;
                    if (fragLen >= minFragLen) {
                        String      name = rec.getName() + "/split" + fragNum + "of" + total;
                        FastqRecord frag = clipFragment(rec, clipFrom + start, fragLen, name);
                        if (frag != null) fragments.add(frag);
                    }
                    fragNum++;
                    start = hit[1] - clipFrom;
                }
                int lastLen = workLen - start;
                if (lastLen >= minFragLen) {
                    String      name = rec.getName() + "/split" + fragNum + "of" + total;
                    FastqRecord frag = clipFragment(rec, clipFrom + start, lastLen, name);
                    if (frag != null) fragments.add(frag);
                }

                if (fragments.isEmpty())   return new FastqRecord[]{null};
                if (fragments.size() == 1) return new FastqRecord[]{fragments.get(0)};
                return fragments.toArray(new FastqRecord[0]);
            }
        }

        if (clipFrom == 0 && trimTo == seqLen) return in;
        return new FastqRecord[]{new FastqRecord(rec, clipFrom, workLen)};
    }

    // -----------------------------------------------------------------------
    // Terminal clipping

    /**
     * Returns {@code seq}'s first {@code len} characters in a reusable per-thread
     * buffer, avoiding a fresh {@code toCharArray()} allocation per call.
     *
     * <p>Each record previously allocated three of these, one in each terminal
     * scan and one in the interior scan, plus one more per emitted fragment via
     * {@link #clipFragment}.  At PacBio HiFi read lengths (15.5 kb average) that is
     * ~31 kB per array, and the class already claims a zero-allocation scan.
     *
     * <p><b>Constraint:</b> the returned array is valid only until the next call on
     * the same thread.  Every current caller consumes it within one method and does
     * not call back into another method that also asks for it, so this is safe; a
     * new caller that needs to hold it across such a call must copy it.
     */
    private char[] readChars(String seq, int len) {
        char[] buf = charScratch.get();
        if (buf.length < len) {
            buf = new char[Math.max(len, buf.length * 2)];
            charScratch.set(buf);
        }
        seq.getChars(0, len, buf, 0);
        return buf;
    }

    /**
     * Edit budget for a candidate match of {@code len} bases.
     *
     * <p>Single source of truth.  The terminal scans, the interior scan and the
     * k-mer pre-filter reliability conditions must all agree on this number: if
     * the pre-filter assumes a smaller budget than the matcher uses, it discards
     * candidates the matcher would have accepted.
     */
    private int allowedEdits(int len) {
        return (int) (len * maxErrorRate);
    }

    /**
     * Returns the number of 5′ bases to remove.
     *
     * <p>Two scans are performed:
     * <ol>
     *   <li><b>Suffix scan</b> (forward adapters only): the adapter may hang off the
     *       read's 5′ end; partial overlaps down to {@code minOverlap} are accepted.</li>
     *   <li><b>Near-terminal full-adapter scan</b> (all orientations): ONT reads
     *       typically have 2–10 bases of pore-entry artifact before the adapter begins.
     *       This scan checks whether a full adapter starts at any position 0..minOverlap,
     *       complementing the interior scan which starts from position minOverlap.</li>
     * </ol>
     */
    private int findFivePrimeClip(String seq, int seqLen) {
        int    clipFrom  = 0;
        char[] readChars = readChars(seq, seqLen);

        // Suffix scan: adapter hangs off the 5′ end (forward adapters only).
        for (int fi = 0; fi < fwdAdapters.size(); fi++) {
            String adapter    = fwdAdapters.get(fi);
            int    adapterLen = adapter.length();
            int    maxOverlap = Math.min(seqLen, adapterLen);
            if (maxOverlap < minOverlap) continue;
            // Pre-filter: safe only when minOverlap >= KMER_SIZE, the adapter has no N
            // wildcards (which leave gaps in the k-mer table), and the overlap is large
            // enough that errors cannot destroy every possible shared k-mer.
            // Reliability condition: (minOverlap - KMER_SIZE + 1) > allowedEdits(minOverlap) * KMER_SIZE
            if (!fwdAdapterHasN[fi] && minOverlap >= KMER_SIZE
                    && (minOverlap - KMER_SIZE + 1) > allowedEdits(minOverlap) * KMER_SIZE
                    && !hasSharedKmer(readChars, 0, maxOverlap, fwdAdapterHasKmer[fi]))
                continue;

            for (int overlap = maxOverlap; overlap >= minOverlap; overlap--) {
                if (overlap <= clipFrom) break;
                int adapterStart = adapterLen - overlap;
                int allowedEdits = allowedEdits(overlap);
                if (editDistance(seq, 0, adapter, adapterStart, overlap, allowedEdits)
                        <= allowedEdits) {
                    clipFrom = overlap;
                    break;
                }
            }
        }

        // Near-terminal full-adapter scan: full adapter starting at positions 0..minOverlap.
        for (int ai = 0; ai < adapters.size(); ai++) {
            String adapter    = adapters.get(ai);
            int    adapterLen = adapter.length();
            int    maxStart   = Math.min(minOverlap, seqLen - adapterLen - minOverlap);
            if (maxStart < 0) continue;
            int filterEnd = Math.min(seqLen, minOverlap + adapterLen);
            int allowedEdits = allowedEdits(adapterLen);
            if (!adapterHasN[ai] && adapterLen >= KMER_SIZE
                    && (adapterLen - KMER_SIZE + 1) > allowedEdits * KMER_SIZE
                    && !hasSharedKmer(readChars, 0, filterEnd, adapterHasKmer[ai]))
                continue;
            for (int s = 0; s <= maxStart; s++) {
                if (editDistance(seq, s, adapter, 0, adapterLen, allowedEdits)
                        <= allowedEdits) {
                    int newClip = s + adapterLen;
                    if (newClip > clipFrom) clipFrom = newClip;
                    break;
                }
            }
        }

        return clipFrom;
    }

    /**
     * Returns the 3′ keep boundary (exclusive): retain seq[0..trimTo-1].
     *
     * <p>One scan is performed: a <b>prefix scan</b> over all orientations, where the
     * adapter may hang off the read's 3′ end and partial overlaps down to
     * {@code minOverlap} are accepted.
     *
     * <p>Adapters that end a few bases before the true 3′ terminus (a trailing barcode
     * or fill bases after the adapter, so the adapter does not hang off the very end)
     * are not this method's responsibility. {@link #findInternalHits} covers them: it
     * runs on all platforms, requires a k-mer seed and the complete adapter, and its
     * hits are reconciled against {@code trimTo} by the caller.
     */
    private int findThreePrimeClip(String seq, int seqLen) {
        int    trimTo    = seqLen;
        char[] readChars = readChars(seq, seqLen);

        // Prefix scan: adapter hangs off the 3′ end.
        for (int ai = 0; ai < adapters.size(); ai++) {
            String adapter    = adapters.get(ai);
            int    adapterLen = adapter.length();
            int    maxOverlap = Math.min(seqLen, adapterLen);
            if (maxOverlap < minOverlap) continue;
            int filterStart = seqLen - maxOverlap;
            // Pre-filter: same safety conditions as the 5′ scan (see findFivePrimeClip).
            if (!adapterHasN[ai] && minOverlap >= KMER_SIZE
                    && (minOverlap - KMER_SIZE + 1) > allowedEdits(minOverlap) * KMER_SIZE
                    && !hasSharedKmer(readChars, filterStart, seqLen, adapterHasKmer[ai]))
                continue;

            for (int overlap = maxOverlap; overlap >= minOverlap; overlap--) {
                int readStart = seqLen - overlap;
                if (readStart >= trimTo) break;
                int allowedEdits = allowedEdits(overlap);
                if (editDistance(seq, readStart, adapter, 0, overlap, allowedEdits)
                        <= allowedEdits) {
                    trimTo = readStart;
                    break;
                }
            }
        }

        // Deliberately no second, near-terminal scan here.  findInternalHits covers
        // that region on every platform, with a k-mer seed and a full-adapter
        // requirement.  A brute-force scan over every start position in a minOverlap
        // window, with no seed requirement, costs precision for nothing: measured on
        // the HiFi ground truth it produced 2486 bases of adjacent over-clipping
        // against 238 without it, at identical recall 1.0000.
        return trimTo;
    }

    /**
     * Re-clips the terminal ends of a split fragment.
     */
    private FastqRecord clipFragment(FastqRecord base, int offset, int len, String name) {
        String fragSeq = base.getSequence().substring(offset, offset + len);
        int    fc      = findFivePrimeClip(fragSeq, len);
        int    ft      = findThreePrimeClip(fragSeq, len);
        int    newLen  = ft - fc;
        if (newLen <= 0) return null;
        return new FastqRecord(base, offset + fc, newLen, name, name);
    }

    // -----------------------------------------------------------------------
    // Interior chimera splitting

    /**
     * Returns sorted, non-overlapping {hitStart, hitEnd} intervals for adapter
     * matches found in the interior of {@code seq}, excluding terminal zones.
     *
     * <p>The k-mer scan uses a flat 4096-entry integer-indexed table
     * ({@link #kmerTable}) and encodes each read position as a 12-bit integer
     * via {@link #encodeKmer}.  This avoids all String allocation in the scan
     * loop. For a 50 kb read that is ~50 000 eliminated String objects per
     * thread compared with the previous {@code seq.substring()} approach.
     *
     * <p>Candidate positions are verified by {@link #bestMatchAt}, which
     * requires the <em>full adapter</em> to match (not just {@code minOverlap}
     * bases).  This prevents false-positive chimera splits: a 10-bp partial
     * interior match at 15 % error tolerates 1 edit and is easily satisfied by
     * chance in high-error-rate reads; requiring the full 20 bp adapter raises
     * the bar to 3 edits in 20 bp, which is far more specific.
     */
    private List<int[]> findInternalHits(String seq, int seqLen) {
        int terminalZone = minOverlap;

        char[] seqChars = readChars(seq, seqLen);
        // Candidates go into a reusable int[] instead of a HashSet<Integer>.  The
        // set allocated one hash table plus one boxed Integer per candidate for
        // every read, which dominated allocation once this scan started running on
        // PacBio HiFi: 3.09M reads at 15.5 kb average.  Sorting then skipping equal
        // neighbours deduplicates, and the sorted order is wanted downstream anyway.
        int[] cand  = candScratch.get();
        int   nCand = 0;
        // Seed to the last position a k-mer fits.  Stopping terminalZone bases
        // early starved the 3' side: an adapter ending d bases from the terminus
        // only kept seeds at adapter offsets <= d + terminalZone - KMER_SIZE, so a
        // near-terminal adapter carrying a couple of errors lost every seed.
        int scanEnd = seqLen - KMER_SIZE;

        for (int rp = terminalZone; rp <= scanEnd; rp++) {
            int code = encodeKmer(seqChars, rp, KMER_SIZE);
            if (code < 0) continue;
            List<int[]> entries = kmerTable[code];
            if (entries == null) continue;
            for (int[] e : entries) {
                int readStart = rp - e[1]; // back-project: seed at adapter position e[1]
                // No 3' upper bound: bestMatchAt already rejects any position where
                // the full adapter does not fit, which is the only real constraint.
                // The old bound (seqLen - terminalZone - minOverlap) excluded every
                // adapter ending within (terminalZone + minOverlap - adapterLen) of
                // the 3' end, a region the terminal prefix scan cannot reach either
                // once the adapter stops hanging off the very last base.
                if (readStart >= terminalZone) {
                    if (nCand == cand.length) {
                        cand = java.util.Arrays.copyOf(cand, cand.length * 2);
                        candScratch.set(cand);
                    }
                    cand[nCand++] = readStart;
                }
            }
        }
        if (nCand == 0) return Collections.emptyList();

        java.util.Arrays.sort(cand, 0, nCand);

        List<int[]> verified = new ArrayList<>();
        int prev = -1;
        for (int i = 0; i < nCand; i++) {
            int readStart = cand[i];
            if (readStart == prev) continue;   // duplicate seeds project to one start
            prev = readStart;
            int[] best = bestMatchAt(seq, seqLen, readStart, terminalZone);
            if (best != null) verified.add(best);
        }
        if (verified.isEmpty()) return Collections.emptyList();

        verified.sort((a, b) -> a[0] - b[0]);
        return mergeOverlapping(verified);
    }

    /**
     * Returns the {hitStart, hitEnd} of the best full-adapter match at
     * {@code readStart}, or {@code null} if no adapter passes the threshold.
     *
     * <p>The minimum acceptable overlap is {@code adapterLen} (not
     * {@code minOverlap}).  Partial interior matches are biologically implausible
     * at a chimera junction and produce false positives in high-error-rate data.
     * If the adapter does not fully fit between {@code readStart} and the scan-zone
     * boundary, the position is skipped.
     */
    private int[] bestMatchAt(String seq, int seqLen, int readStart, int terminalZone) {
        int[] best = null;
        for (String adapter : adapters) {
            int adapterLen = adapter.length();
            // Only accept the full adapter; partial interior matches are rejected.
            // The check is a BOUNDS check only.  It previously also subtracted
            // terminalZone, which excluded any adapter ending within minOverlap
            // of the 3' end, a region the 3' terminal clip cannot reach either,
            // because the surviving overlap there is below minOverlap.
            if (seqLen - readStart < adapterLen) continue;

            int allowedEdits = allowedEdits(adapterLen);
            if (editDistance(seq, readStart, adapter, 0, adapterLen, allowedEdits)
                    <= allowedEdits) {
                if (best == null || adapterLen > (best[1] - best[0]))
                    best = new int[]{readStart, readStart + adapterLen};
            }
        }
        return best;
    }

    /** Merges overlapping or adjacent intervals, preserving sorted order. */
    private List<int[]> mergeOverlapping(List<int[]> sorted) {
        List<int[]> merged  = new ArrayList<>();
        int[]       current = sorted.get(0).clone();
        for (int i = 1; i < sorted.size(); i++) {
            int[] next = sorted.get(i);
            if (next[0] < current[1]) current[1] = Math.max(current[1], next[1]);
            else { merged.add(current); current = next.clone(); }
        }
        merged.add(current);
        return merged;
    }
}
