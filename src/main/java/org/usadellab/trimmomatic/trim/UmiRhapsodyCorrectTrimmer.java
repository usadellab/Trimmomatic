package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * UMIRHAPSODYCORRECT:[<beadVersion>:]<whitelistDir>:<maxMismatch>[:<separator>]
 *
 * Whitelist-corrects the three cell label segments (CLS1/CLS2/CLS3) of a BD
 * Rhapsody R1 read and extracts the combined cell barcode plus the raw UMI to
 * the read name, matching the same bare-sequence convention as
 * UMISPLIT/UMIDROPLETCORRECT. One trimmer covers every supported bead version --
 * <beadVersion> selects the segment lengths, linker lengths and prefix-inset
 * handling internally, the same way LONGREADTRIM's platform hint (ONT/CLR/HIFI)
 * selects behaviour within one class. Neither needs one class per platform.
 *
 * <beadVersion> is optional and defaults to "V1" if omitted. The constructor
 * explains why it is a leading token, unlike LONGREADTRIM's trailing platform
 * hint.
 *
 * ===========================================================================
 * V1 (default; empirically confirmed against real data and BD's own public
 * codebook)
 * ===========================================================================
 * R1 structure, 0-indexed:
 *   CLS1(0-9) - L1(9-21, fixed "ACTGGCCTGCGA") - CLS2(21-30) -
 *   L2(30-43, fixed "GGTAGCGGTGACA") - CLS3(43-52) - UMI(52-60) - poly-T carryover
 * Each CLS is one of 96 known sequences. CLS1 is always read at its fixed
 * nominal position: nothing precedes it to derive an offset from. A single
 * cumulative offset is then derived from CLS2 (checked at nominal, then +/-1,
 * then +/-2, this catches any indel that occurred anywhere before CLS2, in
 * CLS1 or L1, since only the cumulative shift at CLS2's position matters) and
 * applied to CLS3 and the UMI too, without re-deriving it a second time at L2.
 * This was verified empirically on real BD Rhapsody excerpt data (see
 * benchmark/scripts/rhapsody_offset_cascade_analysis.py, not shipped):
 * re-deriving the offset at L2 resolved fewer than 1% of the cases where
 * CLS2's offset failed to place L2 exactly, and moved the CLS3 whitelist-hit
 * rate by under 0.2 percentage points across all three real samples tested.
 * The other ~99% of L2 mismatches are simple substitution noise on a 13bp
 * window, not a second indel, an offset search cannot fix a substitution.
 * If CLS2 has no exact match at any of the five candidate offsets, one
 * Hamming-1 correction attempt is made at the nominal (zero) offset only.
 *
 * ===========================================================================
 * ENHANCEDV2 (best-guess design, NOT validated against real data, see below)
 * ===========================================================================
 * Named after BD's own "Enhanced V2" bead designation specifically, BD's
 * bead line-up isn't a clean V1/V2/V3 progression: "V1" and "Enhanced" are
 * separate bead families, and "Enhanced" itself has sub-generations
 * ("Enhanced", "Enhanced V2", "Enhanced V3"). The structure below is what's
 * documented for "Enhanced V2" specifically, not the whole Enhanced family.
 *
 * R1 structure per BD's own documentation and the community-maintained
 * scg_lib_structs reference (https://teichlab.github.io/scg_lib_structs/methods_html/BD_Rhapsody.html),
 * neither of which this project has cross-checked against actual Enhanced V2
 * reads the way V1 was checked:
 *   [prefix inset: 0-3bp, one of "", "A", "GT", "TCA"] - CLS1(9bp) -
 *   L1(~4bp, "GTGA") - CLS2(9bp) - L2(~4bp, "GACA") - CLS3(9bp) - UMI(8bp) -
 *   poly-T carryover
 * Each CLS is one of 384 known sequences (vs. V1's 96). The prefix inset is a
 * deliberate, always-present variable-length element (not sequencing noise),
 * so it is resolved first. CLS1 is searched for at each candidate inset length
 * (0, 1, 2, 3, in that order, first exact whitelist hit wins) instead of at a
 * single fixed position. Whatever inset length is found becomes the base
 * offset, and CLS2/CLS3/UMI localisation then proceeds exactly as in V1
 * (CLS2 offset search, then straight-line propagation to CLS3/UMI), just
 * shifted by that base offset. If no inset length gives CLS1 an exact match,
 * one Hamming-1 attempt is made at inset length 0.
 *
 * This mode is a documentation-derived best guess, not an empirically
 * verified design: the exact linker sequences/lengths, whether the inset is
 * really limited to those four values, and whether V1's "one offset from CLS2
 * propagates cleanly to CLS3" finding even holds for Enhanced V2 beads are all
 * unconfirmed against real reads. Treat ENHANCEDV2 as experimental until
 * checked against an actual Enhanced V2 dataset.
 *
 * ===========================================================================
 * Common to both versions
 * ===========================================================================
 * The CLS1/CLS2/CLS3 codebooks are BD Biosciences' own commercial-kit
 * reference data with no confirmed redistribution license, so they are not
 * bundled with Trimmomatic, <whitelistDir> must be supplied and must
 * contain CLS1.txt, CLS2.txt and CLS3.txt (one sequence per line, plain text;
 * 96 entries each for V1, 384 each for Enhanced V2, small enough in either
 * case that compressed-format support would add complexity for no gain). A single
 * directory argument is used instead of three separate file paths because
 * three colon-delimited Windows paths in one argument cannot be parsed back
 * apart from each other.
 *
 * A read whose CLS1 (inset search), CLS2 (offset) or CLS3 has no confident
 * match within maxMismatch is dropped, as is a read too short for the
 * selected version's full construct.
 *
 * The read is renamed as:
 *   @original_name<separator><CLS1+CLS2+CLS3><separator><raw UMI>
 *
 * Examples:
 *   UMIRHAPSODYCORRECT:benchmark/barcode_whitelists/rhapsody_v1:1
 *   UMIRHAPSODYCORRECT:ENHANCEDV2:benchmark/barcode_whitelists/rhapsody_enhancedv2:1
 */
public class UmiRhapsodyCorrectTrimmer extends AbstractSingleRecordTrimmer {
    private static final int MAX_SUPPORTED_MISMATCH = 1;
    private static final char[] BASES = { 'A', 'C', 'G', 'T' };

    private static final int[] SEARCH_OFFSETS = { 0, -1, 1, -2, 2 };

    private static final class BeadProfile {
        final int cls1Len, l1Len, cls2Len, l2Len, cls3Len, umiLen;
        final int[] insetCandidates;

        BeadProfile(int cls1Len, int l1Len, int cls2Len, int l2Len, int cls3Len, int umiLen, int[] insetCandidates) {
            this.cls1Len = cls1Len;
            this.l1Len = l1Len;
            this.cls2Len = cls2Len;
            this.l2Len = l2Len;
            this.cls3Len = cls3Len;
            this.umiLen = umiLen;
            this.insetCandidates = insetCandidates;
        }
    }

    // V1: no prefix inset, candidate list of just {0} makes the inset-search
    // loop below degenerate into "always try position 0", identical to this
    // trimmer's original V1-only behaviour.
    private static final BeadProfile PROFILE_V1 = new BeadProfile(9, 12, 9, 13, 9, 8, new int[] { 0 });
    // ENHANCEDV2: best-guess, see class javadoc, NOT validated against real data.
    private static final BeadProfile PROFILE_ENHANCEDV2 = new BeadProfile(9, 4, 9, 4, 9, 8, new int[] { 0, 1, 2, 3 });

    private final int cls1Len, cls2Len, cls3Len, umiLen;
    private final int cls2NominalStart, cls3NominalStart, umiNominalStart, umiNominalEnd;
    private final int[] insetCandidates;

    private final int maxMismatch;
    private final String separator;
    private final Set<String> cls1Whitelist;
    private final Set<String> cls2Whitelist;
    private final Set<String> cls3Whitelist;

    public UmiRhapsodyCorrectTrimmer(String args) throws IOException {
        String[] tokens = args.split(":");
        if (tokens.length < 1)
            throw new IllegalArgumentException(
                    "UMIRHAPSODYCORRECT requires <whitelistDir>:<maxMismatch>, got: " + args);

        // beadVersion is a leading token, not trailing like LONGREADTRIM's
        // platform hint: <separator> is already an optional trailing
        // non-numeric token, and a second optional trailing non-numeric token
        // (beadVersion) couldn't be told apart from it by the
        // is-the-last-token-an-integer heuristic used below. Recognising a
        // fixed vocabulary ("V1"/"ENHANCEDV2") at the front avoids that clash
        // and keeps every existing V1 invocation working unchanged.
        BeadProfile profile;
        String[] remaining;
        if (tokens[0].equalsIgnoreCase("V1")) {
            profile = PROFILE_V1;
            remaining = Arrays.copyOfRange(tokens, 1, tokens.length);
        } else if (tokens[0].equalsIgnoreCase("ENHANCEDV2")) {
            profile = PROFILE_ENHANCEDV2;
            remaining = Arrays.copyOfRange(tokens, 1, tokens.length);
        } else {
            profile = PROFILE_V1;
            remaining = tokens;
        }

        if (remaining.length < 2)
            throw new IllegalArgumentException(
                    "UMIRHAPSODYCORRECT requires <whitelistDir>:<maxMismatch>, got: " + args);

        boolean hasSeparator = !isInteger(remaining[remaining.length - 1]);
        int trailingCount = hasSeparator ? 2 : 1;
        int pathTokenCount = remaining.length - trailingCount;
        if (pathTokenCount < 1)
            throw new IllegalArgumentException(
                    "UMIRHAPSODYCORRECT requires <whitelistDir>:<maxMismatch>, got: " + args);

        String whitelistDir = String.join(":", Arrays.copyOfRange(remaining, 0, pathTokenCount));
        maxMismatch = Integer.parseInt(remaining[pathTokenCount]);
        separator = hasSeparator ? remaining[pathTokenCount + 1] : "_";

        if (maxMismatch < 0 || maxMismatch > MAX_SUPPORTED_MISMATCH)
            throw new IllegalArgumentException(
                    "UMIRHAPSODYCORRECT maxMismatch must be between 0 and " + MAX_SUPPORTED_MISMATCH
                    + " (got " + maxMismatch + "). Larger neighbourhoods cost too much for a plain Hamming correction");

        cls1Len = profile.cls1Len;
        cls2Len = profile.cls2Len;
        cls3Len = profile.cls3Len;
        umiLen = profile.umiLen;
        insetCandidates = profile.insetCandidates;

        cls2NominalStart = cls1Len + profile.l1Len;
        cls3NominalStart = cls2NominalStart + cls2Len + profile.l2Len;
        umiNominalStart = cls3NominalStart + cls3Len;
        umiNominalEnd = umiNominalStart + umiLen;

        File dir = new File(whitelistDir);
        cls1Whitelist = loadWhitelist(new File(dir, "CLS1.txt"));
        cls2Whitelist = loadWhitelist(new File(dir, "CLS2.txt"));
        cls3Whitelist = loadWhitelist(new File(dir, "CLS3.txt"));
    }

    private static boolean isInteger(String s) {
        try {
            Integer.parseInt(s);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private static Set<String> loadWhitelist(File f) throws IOException {
        Set<String> set = new HashSet<>();

        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty())
                    set.add(line);
            }
        }

        if (set.isEmpty())
            throw new IllegalArgumentException("UMIRHAPSODYCORRECT whitelist file is empty or unreadable: " + f);

        return set;
    }

    /**
     * Same symmetric-mode hazard as UMIDROPLETCORRECT, see its processRecords()
     * for the full explanation. Refusing keeps mate names in sync.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMIRHAPSODYCORRECT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len < umiNominalEnd)
            return null; // shorter than the full construct even with zero prefix inset

        String seq = in.getSequence();

        // Resolve the prefix inset (V1: always just candidate 0, degenerating
        // to "CLS1 at its fixed nominal position", identical to V1's original
        // behaviour) and correct CLS1 at the same time.
        Integer insetOffset = null;
        String cls1 = null;
        for (int candidate : insetCandidates) {
            int start = candidate;
            int end = candidate + cls1Len;
            if (start < 0 || end > len)
                continue;
            String window = seq.substring(start, end);
            if (cls1Whitelist.contains(window)) {
                insetOffset = candidate;
                cls1 = window;
                break;
            }
        }
        if (insetOffset == null) {
            int fallback = insetCandidates[0];
            int end = fallback + cls1Len;
            if (end <= len) {
                String corrected = correct(seq.substring(fallback, end), cls1Whitelist);
                if (corrected != null) {
                    insetOffset = fallback;
                    cls1 = corrected;
                }
            }
        }
        if (insetOffset == null)
            return null;

        // CLS2 offset search, shifted by whatever prefix inset was resolved
        // above, otherwise identical to V1's original CLS2-anchored search.
        Integer midOffset = null;
        String cls2 = null;
        for (int delta : SEARCH_OFFSETS) {
            int start = insetOffset + cls2NominalStart + delta;
            int end = start + cls2Len;
            if (start < 0 || end > len)
                continue;
            String window = seq.substring(start, end);
            if (cls2Whitelist.contains(window)) {
                midOffset = delta;
                cls2 = window;
                break;
            }
        }
        if (midOffset == null) {
            int start = insetOffset + cls2NominalStart;
            int end = start + cls2Len;
            if (end <= len) {
                String corrected = correct(seq.substring(start, end), cls2Whitelist);
                if (corrected != null) {
                    midOffset = 0;
                    cls2 = corrected;
                }
            }
        }
        if (midOffset == null)
            return null;

        int totalOffset = insetOffset + midOffset;
        int cls3Start = totalOffset + cls3NominalStart;
        int cls3End = cls3Start + cls3Len;
        int umiStart = totalOffset + umiNominalStart;
        int umiEnd = umiStart + umiLen;
        if (cls3Start < 0 || umiEnd > len)
            return null;

        String cls3 = correct(seq.substring(cls3Start, cls3End), cls3Whitelist);
        if (cls3 == null)
            return null;

        String umi = seq.substring(umiStart, umiEnd);

        // Bare sequences, no "CLS:"/"UMI:" labels, see class javadoc.
        String newName = in.getName() + separator + cls1 + cls2 + cls3 + separator + umi;
        return new FastqRecord(in, umiEnd, len - umiEnd, newName);
    }

    private String correct(String raw, Set<String> whitelist) {
        if (whitelist.contains(raw))
            return raw;

        if (maxMismatch < 1)
            return null;

        char[] chars = raw.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            char original = chars[i];
            for (char base : BASES) {
                if (base == original)
                    continue;
                chars[i] = base;
                String candidate = new String(chars);
                if (whitelist.contains(candidate))
                    return candidate;
            }
            chars[i] = original;
        }

        return null;
    }
}
