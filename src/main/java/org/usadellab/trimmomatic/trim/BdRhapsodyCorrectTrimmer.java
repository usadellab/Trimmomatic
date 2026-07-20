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
 * BDRHAPSODYCORRECT:<whitelistDir>:<maxMismatch>[:<separator>]
 *
 * Whitelist-corrects the three cell label segments (CLS1/CLS2/CLS3) of a BD
 * Rhapsody V1-bead R1 read and extracts the combined 27bp cell barcode plus
 * the raw 8bp UMI to the read name, matching the same bare-sequence
 * convention as UMISPLIT/BARCODECORRECT.
 *
 * R1 structure (empirically confirmed against real data and BD's own public
 * codebook, positions 0-indexed):
 *   CLS1(0-9) - L1(9-21, fixed "ACTGGCCTGCGA") - CLS2(21-30) -
 *   L2(30-43, fixed "GGTAGCGGTGACA") - CLS3(43-52) - UMI(52-60) - poly-T carryover
 * Each CLS is one of 96 known sequences (whitelist-correctable); L1/L2 are
 * fixed linkers used only to anchor position, never extracted; UMI is random
 * and never corrected.
 *
 * Offset handling: a single cumulative offset is derived from CLS2 (checked
 * at nominal, then +/-1, then +/-2 -- this catches any indel that occurred
 * anywhere before CLS2, in CLS1 or L1, since only the cumulative shift at
 * CLS2's position matters). That same offset is then applied to CLS3 and the
 * UMI, without re-deriving it at L2. This was verified empirically on the
 * real BD Rhapsody excerpt data (see benchmark/scripts/rhapsody_offset_cascade_analysis.py,
 * not shipped): re-deriving the offset a second time at L2 resolved fewer
 * than 1% of the cases where CLS2's offset failed to place L2 exactly, and
 * moved the CLS3 whitelist-hit rate by under 0.2 percentage points across all
 * three real samples tested. The other ~99% of L2 mismatches are simple
 * substitution noise on a 13bp window (consistent with typical per-base
 * error rates), not a second indel -- an offset search cannot fix a
 * substitution, so re-checking doesn't help and isn't worth the extra
 * per-read cost.
 *
 * If CLS2 has no exact match at any of the five candidate offsets, one
 * Hamming-1 correction attempt is made at the nominal (zero) offset only --
 * a genuine second indel between CLS1 and CLS2 is rare enough (see above)
 * that a mismatch at every offset is far more likely to be substitution
 * noise at the nominal position than an indel a wider search would catch.
 *
 * CLS1 is always read at its fixed nominal position (0-9): no prior segment
 * exists to derive an offset from before it.
 *
 * The CLS1/CLS2/CLS3 codebooks are BD Biosciences' own commercial-kit
 * reference data (BD Rhapsody V1 bead codebook) with no confirmed
 * redistribution license, so they are not bundled with Trimmomatic --
 * <whitelistDir> must be supplied and must contain CLS1.txt, CLS2.txt and
 * CLS3.txt (one sequence per line, plain text -- 96 tiny entries each for BD
 * Rhapsody V1, unlike the multi-million-line droplet whitelists BARCODECORRECT
 * handles, so compressed-format support isn't worth the complexity here). A
 * single directory argument is used instead of
 * three separate file paths (as BARCODECORRECT takes for its one whitelist)
 * because three colon-delimited Windows paths in one argument cannot be
 * parsed back apart from each other.
 *
 * A read whose CLS1, CLS2 (offset) or CLS3 has no confident match within
 * maxMismatch is dropped. A read shorter than 60bp (the full CLS1-L1-CLS2-L2-
 * CLS3-UMI construct) is dropped outright.
 *
 * The read is renamed as:
 *   @original_name<separator><CLS1+CLS2+CLS3, 27bp><separator><raw UMI, 8bp>
 *
 * Example:
 *   BDRHAPSODYCORRECT:benchmark/barcode_whitelists/rhapsody_v1:1
 */
public class BdRhapsodyCorrectTrimmer extends AbstractSingleRecordTrimmer {
    private static final int MAX_SUPPORTED_MISMATCH = 1;
    private static final char[] BASES = { 'A', 'C', 'G', 'T' };

    private static final int[] SEARCH_OFFSETS = { 0, -1, 1, -2, 2 };

    private static final int CLS1_LEN = 9;
    private static final int L1_LEN = 12;
    private static final int CLS2_LEN = 9;
    private static final int L2_LEN = 13;
    private static final int CLS3_LEN = 9;
    private static final int UMI_LEN = 8;

    private static final int CLS1_START = 0;
    private static final int CLS1_END = CLS1_START + CLS1_LEN;
    private static final int CLS2_START = CLS1_END + L1_LEN;
    private static final int CLS2_END = CLS2_START + CLS2_LEN;
    private static final int CLS3_START = CLS2_END + L2_LEN;
    private static final int CLS3_END = CLS3_START + CLS3_LEN;
    private static final int UMI_START = CLS3_END;
    private static final int UMI_END = UMI_START + UMI_LEN;

    private final int maxMismatch;
    private final String separator;
    private final Set<String> cls1Whitelist;
    private final Set<String> cls2Whitelist;
    private final Set<String> cls3Whitelist;

    public BdRhapsodyCorrectTrimmer(String args) throws IOException {
        // Same right-to-left parsing as BARCODECORRECT (a Windows drive-letter
        // path contains its own colon), but here there is only one path
        // argument -- see class javadoc for why three separate whitelist paths
        // aren't used.
        String[] tokens = args.split(":");
        if (tokens.length < 2)
            throw new IllegalArgumentException(
                    "BDRHAPSODYCORRECT requires <whitelistDir>:<maxMismatch>, got: " + args);

        boolean hasSeparator = !isInteger(tokens[tokens.length - 1]);
        int trailingCount = hasSeparator ? 2 : 1;
        int pathTokenCount = tokens.length - trailingCount;
        if (pathTokenCount < 1)
            throw new IllegalArgumentException(
                    "BDRHAPSODYCORRECT requires <whitelistDir>:<maxMismatch>, got: " + args);

        String whitelistDir = String.join(":", Arrays.copyOfRange(tokens, 0, pathTokenCount));
        maxMismatch = Integer.parseInt(tokens[pathTokenCount]);
        separator = hasSeparator ? tokens[pathTokenCount + 1] : "_";

        if (maxMismatch < 0 || maxMismatch > MAX_SUPPORTED_MISMATCH)
            throw new IllegalArgumentException(
                    "BDRHAPSODYCORRECT maxMismatch must be between 0 and " + MAX_SUPPORTED_MISMATCH
                    + " (got " + maxMismatch + "); larger neighbourhoods aren't worth it for a plain Hamming correction");

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
            throw new IllegalArgumentException("BDRHAPSODYCORRECT whitelist file is empty or unreadable: " + f);

        return set;
    }

    /**
     * Same symmetric-mode hazard as BARCODECORRECT -- see its processRecords()
     * for the full explanation. Refuse rather than desync mate names.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "BDRHAPSODYCORRECT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len < UMI_END)
            return null; // shorter than the full CLS1-L1-CLS2-L2-CLS3-UMI construct

        String seq = in.getSequence();

        // CLS1 is always at the fixed nominal position: nothing precedes it to
        // derive an offset from.
        String cls1 = correct(seq.substring(CLS1_START, CLS1_END), cls1Whitelist);
        if (cls1 == null)
            return null;

        Integer offset = null;
        String cls2 = null;
        for (int delta : SEARCH_OFFSETS) {
            int start = CLS2_START + delta;
            int end = CLS2_END + delta;
            if (start < 0 || end > len)
                continue;
            String window = seq.substring(start, end);
            if (cls2Whitelist.contains(window)) {
                offset = delta;
                cls2 = window;
                break;
            }
        }
        if (offset == null) {
            // No exact match at any offset -- most likely substitution noise at
            // the nominal position rather than an indel a wider search missed
            // (see class javadoc). One Hamming-1 attempt at offset 0 only.
            String corrected = correct(seq.substring(CLS2_START, CLS2_END), cls2Whitelist);
            if (corrected != null) {
                offset = 0;
                cls2 = corrected;
            }
        }
        if (offset == null)
            return null;

        int cls3Start = CLS3_START + offset;
        int cls3End = CLS3_END + offset;
        int umiStart = UMI_START + offset;
        int umiEnd = UMI_END + offset;
        if (cls3Start < 0 || umiEnd > len)
            return null;

        String cls3 = correct(seq.substring(cls3Start, cls3End), cls3Whitelist);
        if (cls3 == null)
            return null;

        String umi = seq.substring(umiStart, umiEnd);

        // Bare sequences, no "CLS:"/"UMI:" labels -- see class javadoc.
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
