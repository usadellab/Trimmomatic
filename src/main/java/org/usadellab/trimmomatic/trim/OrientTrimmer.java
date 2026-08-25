package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * ORIENT:&lt;fasta&gt;:&lt;maxErrorRate&gt;[:&lt;onMissing&gt;]
 *
 * <p>Puts long reads on a common strand. It locates a known fixed sequence and
 * reverse-complements the reads that carry it on the opposite strand.
 *
 * <p>Nanopore and PacBio cDNA libraries are sequenced in both orientations.
 * About half the reads carry their 5' structure reverse-complemented at the 3'
 * end. Any step that reads structure from a fixed position, including most of
 * the UMI family, finds nothing at the expected offset in those reads and drops
 * them. On a scCOLOR-seq run (SRR28589563) the library's fixed anchor appears
 * forward in 64.5% of reads and reverse-complemented in 33.9%. Without this
 * step a third of the data is lost.
 *
 * <p>Orientation is its own task, separate from clipping and extraction, so it
 * is its own step. Put it before the steps that need a common strand and they
 * work unchanged:
 * <pre>
 *   ORIENT:anchor.fa:0.15  LONGREADTRIM:adapters.fa:0.10  HEADCROP:2  UMIDIMERCORRECT:whitelist.txt:24:16:1
 * </pre>
 *
 * <p><b>Matching.</b> Each FASTA entry is searched in the read in both
 * orientations. Candidate positions come from shared 6-mers. Each candidate is
 * confirmed with the same bit-parallel edit distance LONGREADTRIM uses, so an
 * indel costs one edit. A match counts when its edit distance is at most
 * {@code floor(queryLength * maxErrorRate)}.
 *
 * <p>The orientation with the lower edit distance wins. A forward win emits the
 * read unchanged. A reverse win emits the read reverse-complemented, with bases
 * complemented and reversed, quality string reversed, name and comment kept.
 *
 * <p><b>Reads that cannot be oriented.</b> This covers two cases. First, no
 * match within budget in either orientation. Second, an exact tie between the
 * two, which happens on chimeric reads that carry the anchor twice. Both follow
 * {@code onMissing}. {@code DROP} (default) discards the read. {@code KEEP}
 * passes it through untouched, for pipelines where a later step copes with
 * mixed orientation.
 *
 * <p><b>Single-end only.</b> Reverse-complementing one mate of a pair without
 * the other breaks the pairing. No current long-read platform produces
 * paired-end data, so paired-end mode is refused.
 *
 * <p>Example:
 * <pre>
 *   ORIENT:smart_primer.fa:0.15
 * </pre>
 * A FASTA holding the SMART primer AAGCAGTGGTATCAACGCAGAGT. This allows
 * {@code floor(23 * 0.15) = 3} edits, the same tolerance the scCOLOR-seq
 * authors' own tooling uses for that primer.
 */
public class OrientTrimmer extends AbstractSingleRecordTrimmer {
    private static final int KMER_SIZE = 6;
    private static final int KMER_TABLE_SIZE = 1 << (2 * KMER_SIZE);
    private static final int MAX_QUERY_LENGTH = 64;

    /** What to do with a read whose orientation cannot be decided. */
    public enum OnMissing { DROP, KEEP }

    private final List<String> fwdQueries = new ArrayList<>();
    private final List<String> revQueries = new ArrayList<>();

    /** 6-mer seed index, one per orientation: code to list of {queryIndex, offsetInQuery}. */
    @SuppressWarnings("unchecked")
    private final List<int[]>[] fwdSeeds = new List[KMER_TABLE_SIZE];
    @SuppressWarnings("unchecked")
    private final List<int[]>[] revSeeds = new List[KMER_TABLE_SIZE];

    private final float maxErrorRate;
    private final OnMissing onMissing;

    public OrientTrimmer(String args) throws IOException {
        String[] parts = args.split(":");
        if (parts.length < 2)
            throw new IllegalArgumentException(
                    "ORIENT requires <fasta>:<maxErrorRate>, got: " + args);

        // The FASTA path can contain colons (a Windows drive letter). Peel the
        // trailing tokens off the end, the same way UMIDROPLETCORRECT and
        // UMIDIMERCORRECT parse their paths.
        int trailing = 1;
        OnMissing tempOnMissing = OnMissing.DROP;
        String last = parts[parts.length - 1].trim().toUpperCase();
        if (last.equals("DROP") || last.equals("KEEP")) {
            tempOnMissing = OnMissing.valueOf(last);
            trailing = 2;
        }
        onMissing = tempOnMissing;

        int pathTokens = parts.length - trailing;
        if (pathTokens < 1)
            throw new IllegalArgumentException(
                    "ORIENT requires <fasta>:<maxErrorRate>, got: " + args);

        String fastaPath = String.join(":", java.util.Arrays.copyOfRange(parts, 0, pathTokens));
        maxErrorRate = Float.parseFloat(parts[pathTokens]);
        if (maxErrorRate < 0.0f || maxErrorRate > 1.0f)
            throw new IllegalArgumentException(
                    "ORIENT maxErrorRate must be between 0.0 and 1.0, got: " + maxErrorRate);

        loadQueries(fastaPath);
        if (fwdQueries.isEmpty())
            throw new IllegalArgumentException("ORIENT FASTA file contains no sequences: " + fastaPath);
    }

    private void loadQueries(String path) throws IOException {
        List<String> seqs = new ArrayList<>();
        StringBuilder current = new StringBuilder();

        try (BufferedReader br = new BufferedReader(new FileReader(new File(path)))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty())
                    continue;
                if (line.charAt(0) == '>') {
                    if (current.length() > 0) {
                        seqs.add(current.toString());
                        current.setLength(0);
                    }
                } else {
                    current.append(line);
                }
            }
        }
        if (current.length() > 0)
            seqs.add(current.toString());

        for (String raw : seqs) {
            String seq = raw.toUpperCase();
            if (seq.length() < KMER_SIZE)
                throw new IllegalArgumentException(
                        "ORIENT query sequences must be at least " + KMER_SIZE + "bp, got: " + seq);
            // The bit-parallel edit distance is exact only up to 64bp. This step
            // wants a short fixed landmark such as a primer. An anchor longer
            // than 64bp is a configuration mistake.
            if (seq.length() > MAX_QUERY_LENGTH)
                throw new IllegalArgumentException(
                        "ORIENT query sequences must be at most " + MAX_QUERY_LENGTH + "bp, got " + seq.length()
                        + "bp. Supply a short fixed anchor such as a primer, not a full cassette");

            int idx = fwdQueries.size();
            fwdQueries.add(seq);
            revQueries.add(LongReadTrimmer.reverseComplement(seq));
            index(fwdSeeds, fwdQueries.get(idx), idx);
            index(revSeeds, revQueries.get(idx), idx);
        }
    }

    private static void index(List<int[]>[] table, String seq, int queryIndex) {
        char[] chars = seq.toCharArray();
        for (int pos = 0; pos <= chars.length - KMER_SIZE; pos++) {
            int code = LongReadTrimmer.encodeKmer(chars, pos, KMER_SIZE);
            if (code < 0)
                continue; // ambiguous base, cannot seed on it
            if (table[code] == null)
                table[code] = new ArrayList<>();
            table[code].add(new int[] { queryIndex, pos });
        }
    }

    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "ORIENT cannot run in paired-end mode: reverse-complementing one mate without the other "
                    + "would break the pair. Route it to one read only, e.g. via -pe1steps/-pe2steps.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        String seq = in.getSequence();

        int fwdBest = bestDistance(seq, fwdQueries, fwdSeeds);
        int revBest = bestDistance(seq, revQueries, revSeeds);

        boolean fwdHit = fwdBest >= 0;
        boolean revHit = revBest >= 0;

        if (!fwdHit && !revHit)
            return onMissing == OnMissing.KEEP ? in : null;
        if (fwdHit && revHit && fwdBest == revBest)
            return onMissing == OnMissing.KEEP ? in : null;

        if (fwdHit && (!revHit || fwdBest < revBest))
            return in;

        return reverseComplementRecord(in, seq);
    }

    /**
     * Best edit distance of any query against the read, or -1 when nothing
     * matched within budget. Candidates come from shared 6-mers. This keeps the
     * expensive verification to a handful of calls per read. Checking every read
     * position would cost hundreds.
     */
    private int bestDistance(String seq, List<String> queries, List<int[]>[] seeds) {
        int readLen = seq.length();
        if (readLen < KMER_SIZE)
            return -1;

        char[] readChars = seq.toCharArray();
        int best = -1;

        // Marks candidate start positions already verified, so a query matching
        // on several of its own k-mers is not re-checked once per seed.
        boolean[] tried = new boolean[readLen];

        for (int rp = 0; rp <= readLen - KMER_SIZE; rp++) {
            int code = LongReadTrimmer.encodeKmer(readChars, rp, KMER_SIZE);
            if (code < 0)
                continue;
            List<int[]> hits = seeds[code];
            if (hits == null)
                continue;

            for (int[] hit : hits) {
                int qi = hit[0];
                int start = rp - hit[1];
                if (start < 0 || tried[start])
                    continue;

                String query = queries.get(qi);
                int qLen = query.length();
                if (start + qLen > readLen)
                    continue;

                tried[start] = true;
                int budget = (int) (qLen * maxErrorRate);
                int dist = LongReadTrimmer.editDistanceBP(query, 0, seq, start, qLen, budget);
                if (dist <= budget && (best < 0 || dist < best))
                    best = dist;
                if (best == 0)
                    return 0; // exact hit, nothing can beat it
            }
        }

        return best;
    }

    private static FastqRecord reverseComplementRecord(FastqRecord in, String seq) {
        String rcSeq = LongReadTrimmer.reverseComplement(seq);
        String quality = in.getQuality();
        StringBuilder rcQual = new StringBuilder(quality.length());
        for (int i = quality.length() - 1; i >= 0; i--)
            rcQual.append(quality.charAt(i));
        return new FastqRecord(in, rcSeq, rcQual.toString(), in.getPhredOffset());
    }
}
