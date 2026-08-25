package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Set;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.compression.CompressionFormat;

/**
 * UMIDIMERCORRECT:<whitelistFile>:<cbLength>:<umiLength>:<maxMismatchDimers>[:<separator>]
 *
 * Whitelist-corrects a cell barcode built from dimer blocks, the scBUC-seq /
 * scCOLOR-seq design (Philpott et al., Nat Biotechnol 2021,
 * https://doi.org/10.1038/s41587-021-00965-w; see also
 * https://doi.org/10.1101/2021.01.18.427145) used to make barcode/UMI
 * assignment robust to Nanopore's high raw error rate for direct single-cell
 * Nanopore transcriptome sequencing. Each barcode position is synthesised as
 * a 2-base dimer symbol drawn from a restricted set (e.g. only WW/SS-class
 * dimers) instead of a single random base. A single-base sequencing error
 * therefore corrupts at most one dimer symbol and leaves the barcode's read
 * frame intact. A dimer that absorbs one or two base errors still round-trips
 * to the intended symbol far more often than a plain single-base whitelist
 * scheme manages.
 *
 * <cbLength> and <umiLength> are given in bases (not dimers) to stay
 * consistent with every other UMI* trimmer's argument convention; both must
 * be even, since every position is a 2-base dimer by construction.
 *
 * Correction: exact whitelist match on the full barcode first. Failing that,
 * with maxMismatchDimers >= 1, every "single dimer replaced" neighbour of the
 * raw barcode is checked against the whitelist (numDimers positions x 16
 * candidate replacement dimers each). This corrects a barcode whose error
 * is confined to one dimer, regardless of whether that dimer absorbed one or
 * two base-level substitutions, which is exactly the failure mode the dimer
 * design targets. maxMismatchDimers above 1 is rejected outright: two
 * independently-corrupted dimers is a combinatorially large neighbourhood
 * (numDimers choose 2 x 256) for a plain nearest-whitelist-entry correction,
 * mirroring UMIDROPLETCORRECT/UMIRHAPSODYCORRECT's own single-mismatch cap,
 * for the same reason. It costs too much for a correction this simple.
 *
 * The UMI itself is never corrected here, same rationale as
 * UMIDROPLETCORRECT: real UMI correction needs reads grouped by (cell, gene)
 * after alignment, which a pre-alignment trimmer does not have.
 *
 * The read is renamed as:
 *   @original_name<separator><corrected CB bases><separator><raw UMI bases>
 * Bare sequences, matching every other UMI*-family trimmer's umi_tools
 * extract convention.
 *
 * A read whose barcode has no confident match within maxMismatchDimers, or
 * that is shorter than cbLength + umiLength, is dropped.
 *
 * INPUT PRECONDITION: the barcode and UMI must be the leading bases of the
 * read. Raw scCOLOR-seq nanopore reads do NOT meet this. They are full-length
 * cDNA in mixed orientation. The barcode sits behind a library adapter, the
 * SMART primer AAGCAGTGGTATCAACGCAGAGT, and a 2bp spacer. On SRR28589563 that
 * primer occurs around position 31 to 40 in most reads, forward in 64.5% and
 * reverse-complemented in 33.9%. This step has no way to check the
 * precondition. It corrects whatever the leading bases contain, so raw reads
 * produce output that looks valid and is wrong. Put ORIENT and LONGREADTRIM in
 * front of it:
 *
 *   ORIENT:smart_primer.fa:0.15 LONGREADTRIM:adapters.fa:0.10 HEADCROP:2 \
 *   UMIDIMERCORRECT:sccolor_whitelist.txt:24:16:1
 *
 * Example:
 *   UMIDIMERCORRECT:sccolor_whitelist.txt:24:16:1
 *     - 12-dimer (24bp) cell barcode, 16bp raw UMI, 1 dimer-level correction.
 *       The 16bp UMI length is what the scCOLOR-seq authors' own tooling reads
 *       (TallyNN identify_perfect_nano.py takes bases 26 to 42 after the primer,
 *       giving a 40bp barcode-plus-UMI record).
 */
public class UmiDimerCorrectTrimmer extends AbstractSingleRecordTrimmer {
    private static final int MAX_SUPPORTED_MISMATCH_DIMERS = 1;
    private static final char[] BASES = { 'A', 'C', 'G', 'T' };

    private final int cbLength;
    private final int umiLength;
    private final int maxMismatchDimers;
    private final String separator;
    private final Set<String> whitelist;

    public UmiDimerCorrectTrimmer(String args) throws IOException {
        // Same right-to-left token parsing as UMIDROPLETCORRECT. The whitelist
        // path can itself contain colons (a Windows drive letter), so the trailing
        // <cbLength>:<umiLength>:<maxMismatchDimers>[:<separator>] tokens are peeled
        // off from the end instead of splitting the whole string naively.
        String[] tokens = args.split(":");
        if (tokens.length < 4)
            throw new IllegalArgumentException(
                    "UMIDIMERCORRECT requires <whitelistFile>:<cbLength>:<umiLength>:<maxMismatchDimers>, got: " + args);

        boolean hasSeparator = !isInteger(tokens[tokens.length - 1]);
        int trailingCount = hasSeparator ? 4 : 3;
        int pathTokenCount = tokens.length - trailingCount;
        if (pathTokenCount < 1)
            throw new IllegalArgumentException(
                    "UMIDIMERCORRECT requires <whitelistFile>:<cbLength>:<umiLength>:<maxMismatchDimers>, got: " + args);

        String whitelistPath = String.join(":", java.util.Arrays.copyOfRange(tokens, 0, pathTokenCount));
        cbLength = Integer.parseInt(tokens[pathTokenCount]);
        umiLength = Integer.parseInt(tokens[pathTokenCount + 1]);
        maxMismatchDimers = Integer.parseInt(tokens[pathTokenCount + 2]);
        separator = hasSeparator ? tokens[pathTokenCount + 3] : "_";

        if (cbLength < 2 || cbLength % 2 != 0)
            throw new IllegalArgumentException("UMIDIMERCORRECT cbLength must be a positive even number of bases (whole dimers), got: " + cbLength);
        if (umiLength < 1)
            throw new IllegalArgumentException("UMIDIMERCORRECT umiLength must be >= 1, got: " + umiLength);
        if (maxMismatchDimers < 0 || maxMismatchDimers > MAX_SUPPORTED_MISMATCH_DIMERS)
            throw new IllegalArgumentException(
                    "UMIDIMERCORRECT maxMismatchDimers must be between 0 and " + MAX_SUPPORTED_MISMATCH_DIMERS
                    + " (got " + maxMismatchDimers + "). A larger neighbourhood costs too much for a plain nearest-whitelist-entry correction");

        whitelist = loadWhitelist(whitelistPath);
    }

    private static boolean isInteger(String s) {
        try {
            Integer.parseInt(s);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    private static Set<String> loadWhitelist(String path) throws IOException {
        File f = new File(path);
        Set<String> set = new HashSet<>();

        try (InputStream raw = new FileInputStream(f);
                InputStream decompressed = CompressionFormat.wrapStreamForParsing(raw, f.getName());
                BufferedReader br = new BufferedReader(new InputStreamReader(decompressed))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty())
                    set.add(line);
            }
        }

        if (set.isEmpty())
            throw new IllegalArgumentException("UMIDIMERCORRECT whitelist file is empty or unreadable: " + path);

        return set;
    }

    /**
     * Same symmetric-mode hazard as UMIDROPLETCORRECT. See its processRecords()
     * for the full explanation. Refusing keeps mate names in sync.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMIDIMERCORRECT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        int totalLength = cbLength + umiLength;
        if (len < totalLength)
            return null;

        String seq = in.getSequence();
        String rawCb = seq.substring(0, cbLength);
        String umi = seq.substring(cbLength, totalLength);

        String correctedCb = correct(rawCb);
        if (correctedCb == null)
            return null;

        String newName = in.getName() + separator + correctedCb + separator + umi;
        return new FastqRecord(in, totalLength, len - totalLength, newName);
    }

    private String correct(String rawCb) {
        if (whitelist.contains(rawCb))
            return rawCb;

        if (maxMismatchDimers < 1)
            return null;

        int numDimers = rawCb.length() / 2;
        char[] chars = rawCb.toCharArray();

        for (int d = 0; d < numDimers; d++) {
            int pos = d * 2;
            char orig0 = chars[pos];
            char orig1 = chars[pos + 1];

            for (char b0 : BASES) {
                for (char b1 : BASES) {
                    if (b0 == orig0 && b1 == orig1)
                        continue;
                    chars[pos] = b0;
                    chars[pos + 1] = b1;
                    String candidate = new String(chars);
                    if (whitelist.contains(candidate)) {
                        return candidate;
                    }
                }
            }
            chars[pos] = orig0;
            chars[pos + 1] = orig1;
        }

        return null;
    }
}
