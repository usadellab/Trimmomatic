package org.usadellab.trimmomatic.trim;

import java.util.HashMap;
import java.util.Map;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * UMILONGREADEXTRACT:<anchor>:<umiPattern>:<maxMismatch>:<maxIndelShift>[:<separator>]
 *
 * Generic UMI extractor for the long-read amplicon UMI cassette design shared
 * by Oxford Nanopore's own "Custom PCR UMI" protocol (SQK-LSK109) and the
 * equivalent PacBio CCS amplicon-UMI workflow (Karst et al. 2021,
 * https://doi.org/10.1038/s41592-020-01041-y). Both platforms sequence the
 * SAME synthetic cassette (gene-specific primer + fixed anchor + structured
 * degenerate UMI block), just with different downstream basecallers and error
 * profiles. One trimmer therefore covers both platforms.
 *
 * This is intentionally NOT a hardcoded per-kit profile (unlike
 * UMIRHAPSODYCORRECT's V1/ENHANCEDV2 bead profiles). The ONT/PacBio UMI
 * cassette is not one fixed layout across labs and kits: anchor sequence,
 * UMI block length and degeneracy pattern all vary by protocol. Every
 * structural element is a parameter instead.
 *
 * <anchor>: fixed sequence immediately preceding the UMI block (the fixed
 * portion of the tagging primer, e.g. ONT's own UVP-adjacent primer stub).
 * Pass the literal token "NONE" to skip anchor search entirely and assume the
 * UMI starts at position 0. This is required for kits like PCS114/PCB114,
 * where the strand-switching-incorporated UMI cassette sequence is not
 * publicly disclosed by ONT, so no anchor can be hardcoded or supplied.
 *
 * <umiPattern>: an IUPAC degeneracy pattern, one character per UMI base,
 * matching the notation ONT and Karst et al. use for their own cassette
 * designs (e.g. "TTVVVVTTVVVVTTVVVVTTVVVVTTT", V=A/C/G, avoiding the T
 * homopolymer runs that dominate nanopore error). Supported symbols: A/C/G/T
 * (literal, must match exactly), N (any base), and the standard IUPAC
 * 2-fold/3-fold ambiguity codes R/Y/S/W/K/M/B/D/H/V. Pattern length fixes the
 * UMI block length.
 *
 * <maxMismatch>: Hamming mismatch budget, applied twice: once when locating
 * the anchor (plain substitution count against the literal anchor sequence),
 * and once when validating the UMI block against <umiPattern> (a base outside
 * a position's allowed set counts as one mismatch). This is one shared budget.
 * Two separate knobs would double the argument list, and no documented case
 * exists where real users need the two values to differ.
 *
 * <maxIndelShift>: nanopore reads are indel-dominated, not
 * substitution-dominated like Illumina, so a fixed-offset anchor search would
 * misfire on nearly every read whenever a small indel occurs upstream of the
 * anchor. Reuses UMIRHAPSODYCORRECT's offset-cascade technique: the anchor is
 * searched at offset 0 first, then +/-1, +/-2, ... up to +/-maxIndelShift
 * (closest-to-nominal first), and the first offset giving a mismatch count
 * within budget wins. This absorbs any single indel that occurred before the
 * anchor without needing full banded alignment.
 *
 * <separator>: defaults to "_". The UMI is appended to the read name as
 * "<separator>UMI:<bases>", matching UMIEXTRACT's own tagging convention
 * (this trimmer is an "extract" family member, not a "split"/"correct" one,
 * since there is no whitelist to correct a random UMI against).
 *
 * Single 5' end only in this version. The real ONT/Karst cassette can place a
 * UMI at both ends via distinct fwd/rev primers. That needs strand orientation
 * detection, which this trimmer does not perform. Put ORIENT in front of it for
 * mixed-orientation data. If a 3' UMI is also present, run a second instance
 * against a reverse-complemented copy of the data.
 *
 * A read whose anchor cannot be located within budget, or whose UMI block
 * fails pattern validation within budget, is dropped.
 *
 * Examples:
 *   UMILONGREADEXTRACT:GTATCGTGTAGAGACTGCGTAGGT:TTVVVVTTVVVVTTVVVVTTVVVVTTT:2:2
 *     - ONT/Karst-style cassette: fixed anchor, 27bp structured UMI, up to 2
 *       mismatches at anchor/UMI, anchor search shifted up to +/-2bp for indels.
 *   UMILONGREADEXTRACT:NONE:NNNNNNNNNNNNNNNNNN:1:0
 *     - no known anchor (e.g. undisclosed PCS114/PCB114 cassette): 18bp raw N
 *       block assumed at position 0, 1 mismatch tolerated against N (i.e. none,
 *       since N accepts everything, effectively just a length check).
 */
public class UmiLongReadExtractTrimmer extends AbstractSingleRecordTrimmer {
    private static final Map<Character, String> IUPAC = new HashMap<>();
    static {
        IUPAC.put('A', "A");
        IUPAC.put('C', "C");
        IUPAC.put('G', "G");
        IUPAC.put('T', "T");
        IUPAC.put('N', "ACGT");
        IUPAC.put('R', "AG");
        IUPAC.put('Y', "CT");
        IUPAC.put('S', "GC");
        IUPAC.put('W', "AT");
        IUPAC.put('K', "GT");
        IUPAC.put('M', "AC");
        IUPAC.put('B', "CGT");
        IUPAC.put('D', "AGT");
        IUPAC.put('H', "ACT");
        IUPAC.put('V', "ACG");
    }

    private final String anchor;
    private final boolean hasAnchor;
    private final String umiPattern;
    private final int umiLength;
    private final int maxMismatch;
    private final int maxIndelShift;
    private final String separator;

    public UmiLongReadExtractTrimmer(String args) {
        String[] arg = args.split(":");
        if (arg.length < 4)
            throw new IllegalArgumentException(
                    "UMILONGREADEXTRACT requires <anchor>:<umiPattern>:<maxMismatch>:<maxIndelShift>, got: " + args);

        hasAnchor = !arg[0].equalsIgnoreCase("NONE");
        anchor = hasAnchor ? arg[0].toUpperCase() : "";
        if (hasAnchor)
            validateBases(anchor, "anchor");

        umiPattern = arg[1].toUpperCase();
        if (umiPattern.isEmpty())
            throw new IllegalArgumentException("UMILONGREADEXTRACT umiPattern must not be empty");
        for (char c : umiPattern.toCharArray())
            if (!IUPAC.containsKey(c))
                throw new IllegalArgumentException("UMILONGREADEXTRACT umiPattern contains unsupported symbol: " + c);
        umiLength = umiPattern.length();

        maxMismatch = Integer.parseInt(arg[2]);
        if (maxMismatch < 0)
            throw new IllegalArgumentException("UMILONGREADEXTRACT maxMismatch must be >= 0, got: " + maxMismatch);

        maxIndelShift = Integer.parseInt(arg[3]);
        if (maxIndelShift < 0)
            throw new IllegalArgumentException("UMILONGREADEXTRACT maxIndelShift must be >= 0, got: " + maxIndelShift);
        if (maxIndelShift > 5)
            throw new IllegalArgumentException(
                    "UMILONGREADEXTRACT maxIndelShift must be <= 5, got: " + maxIndelShift
                    + ". A wider shift search costs too much, since every offset is re-scanned per read");

        separator = arg.length > 4 ? arg[4] : "_";
    }

    private static void validateBases(String s, String label) {
        for (char c : s.toCharArray())
            if ("ACGT".indexOf(c) < 0)
                throw new IllegalArgumentException("UMILONGREADEXTRACT " + label + " must be A/C/G/T only, got: " + s);
    }

    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMILONGREADEXTRACT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        String seq = in.getSequence();
        int len = seq.length();

        int anchorLen = hasAnchor ? anchor.length() : 0;

        Integer umiStart = null;
        if (!hasAnchor) {
            umiStart = 0;
        } else {
            // Walk 0, -1, +1, -2, +2, ... up to +/-maxIndelShift, closest-to-nominal
            // first, absorbing any single indel that occurred before the anchor.
            for (int i = 0; i <= 2 * maxIndelShift && umiStart == null; i++) {
                int delta = (i == 0) ? 0 : ((i % 2 == 1) ? -((i + 1) / 2) : (i / 2));
                int start = delta;
                int end = start + anchorLen;
                if (start < 0 || end > len)
                    continue;
                if (hamming(seq, start, anchor) <= maxMismatch)
                    umiStart = end;
            }
        }

        if (umiStart == null)
            return null;

        int umiEnd = umiStart + umiLength;
        if (umiEnd > len)
            return null;

        String umi = seq.substring(umiStart, umiEnd);
        if (patternMismatches(umi) > maxMismatch)
            return null;

        String newName = in.getName() + separator + "UMI:" + umi;
        return new FastqRecord(in, umiEnd, len - umiEnd, newName);
    }

    private static int hamming(String seq, int offset, String pattern) {
        int mism = 0;
        for (int i = 0; i < pattern.length(); i++)
            if (seq.charAt(offset + i) != pattern.charAt(i))
                mism++;
        return mism;
    }

    private int patternMismatches(String umi) {
        int mism = 0;
        for (int i = 0; i < umiLength; i++) {
            String allowed = IUPAC.get(umiPattern.charAt(i));
            if (allowed.indexOf(umi.charAt(i)) < 0)
                mism++;
        }
        return mism;
    }
}
