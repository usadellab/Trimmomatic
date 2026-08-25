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
 * UMIDROPLETCORRECT:<whitelistFile>:<cbLength>:<umiLength>:<maxMismatch>[:<separator>]
 *
 * Whitelist-corrects a cell barcode of <cbLength> bases against a known-good
 * barcode list (one sequence per line, plain text or .gz/.bz2/.zip), then
 * appends the CORRECTED barcode and the raw (uncorrected) UMI of
 * <umiLength> bases to the read name, trimming both off the sequence.
 * Combines whitelist correction and header extraction in one step for
 * droplet single-cell platforms (10x Genomics and similar) where the two
 * always happen together.
 *
 * Correction: exact match against the whitelist first; if that fails and
 * maxMismatch >= 1, every 1-substitution neighbour of the raw barcode
 * (cbLength positions x 3 alternate bases, or all 4 at an N) is checked
 * against the whitelist, first hit wins. This is a plain Hamming-distance
 * correction, not the quality- and abundance-weighted Bayesian posterior
 * Cell Ranger's own algorithm uses, it will land in the same ballpark, not
 * be bit-identical to Cell Ranger's numbers. maxMismatch above 1 is rejected
 * outright. The neighbourhood size grows combinatorially, reaching 1080
 * candidates at 2 mismatches for a 16bp barcode, which costs too much for a
 * correction this simple.
 *
 * A read whose barcode has no confident match within maxMismatch is
 * dropped. No leftover payload is required, a read that is exactly
 * cbLength + umiLength bases (the common case: 10x Chromium R1 is 100%
 * barcode+UMI, zero biological payload by design) survives with a
 * zero-length sequence; it is dropped only if shorter than that.
 *
 * The read is renamed as:
 *   @original_name<separator><corrected CB bases><separator><raw UMI bases>
 * These are bare sequences, no "CB:"/"UMI:" labels, matching umi_tools extract's own
 * convention so umi_tools dedup/count's default read-name parser (last
 * underscore-delimited, expects a bare nucleotide string) can consume this
 * directly.
 *
 * UMI is never corrected here, real UMI correction needs reads grouped by
 * (cell, gene) after alignment, which a pre-alignment trimmer does not have.
 *
 * Example:
 *   UMIDROPLETCORRECT:3M-february-2018.txt.gz:16:12:1
 *     - 10x v3 chemistry: 16bp CB (whitelist-corrected, <=1 mismatch), 12bp raw UMI
 */
public class UmiDropletCorrectTrimmer extends AbstractSingleRecordTrimmer {
    private static final int MAX_SUPPORTED_MISMATCH = 1;
    private static final char[] BASES = { 'A', 'C', 'G', 'T' };

    private final int cbLength;
    private final int umiLength;
    private final int maxMismatch;
    private final String separator;
    private final Set<String> whitelist;

    public UmiDropletCorrectTrimmer(String args) throws IOException {
        // The whitelist path can itself contain colons (a Windows drive letter,
        // "C:\...", is the common real case), which breaks a naive split(":")
        // into <path>:<cbLength>:<umiLength>:<maxMismatch>[:<separator>], the
        // path would get sliced apart along with the real args. Parse from the
        // right instead: the last 3 tokens are always cbLength:umiLength:
        // maxMismatch; a 4th trailing token, if present, is the separator,
        // distinguished from "the path had another colon" by whether the last
        // token parses as an integer (a separator string realistically never
        // will, this heuristic would misfire only for a purely numeric
        // separator, which no existing UMIEXTRACT/UMISPLIT test or example uses).
        String[] tokens = args.split(":");
        if (tokens.length < 4)
            throw new IllegalArgumentException(
                    "UMIDROPLETCORRECT requires <whitelistFile>:<cbLength>:<umiLength>:<maxMismatch>, got: " + args);

        boolean hasSeparator = !isInteger(tokens[tokens.length - 1]);
        int trailingCount = hasSeparator ? 4 : 3;
        int pathTokenCount = tokens.length - trailingCount;
        if (pathTokenCount < 1)
            throw new IllegalArgumentException(
                    "UMIDROPLETCORRECT requires <whitelistFile>:<cbLength>:<umiLength>:<maxMismatch>, got: " + args);

        String whitelistPath = String.join(":", java.util.Arrays.copyOfRange(tokens, 0, pathTokenCount));
        cbLength = Integer.parseInt(tokens[pathTokenCount]);
        umiLength = Integer.parseInt(tokens[pathTokenCount + 1]);
        maxMismatch = Integer.parseInt(tokens[pathTokenCount + 2]);
        separator = hasSeparator ? tokens[pathTokenCount + 3] : "_";

        if (cbLength < 1)
            throw new IllegalArgumentException("UMIDROPLETCORRECT cbLength must be >= 1, got: " + cbLength);
        if (umiLength < 1)
            throw new IllegalArgumentException("UMIDROPLETCORRECT umiLength must be >= 1, got: " + umiLength);
        if (maxMismatch < 0 || maxMismatch > MAX_SUPPORTED_MISMATCH)
            throw new IllegalArgumentException(
                    "UMIDROPLETCORRECT maxMismatch must be between 0 and " + MAX_SUPPORTED_MISMATCH
                    + " (got " + maxMismatch + "). Larger neighbourhoods cost too much for a plain Hamming correction");

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
            throw new IllegalArgumentException("UMIDROPLETCORRECT whitelist file is empty or unreadable: " + path);

        return set;
    }

    /**
     * Same symmetric-mode hazard as UMIEXTRACT, see its processRecords() for
     * the full explanation. Refusing keeps mate names in sync.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMIDROPLETCORRECT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        int totalLength = cbLength + umiLength;
        if (len < totalLength)
            return null; // shorter than the barcode+UMI construct itself

        String seq = in.getSequence();
        String rawCb = seq.substring(0, cbLength);
        String umi = seq.substring(cbLength, totalLength);

        String correctedCb = correct(rawCb);
        if (correctedCb == null)
            return null; // no confident whitelist match within maxMismatch

        // Bare sequences, no "CB:"/"UMI:" labels, see class javadoc.
        String newName = in.getName() + separator + correctedCb + separator + umi;
        return new FastqRecord(in, totalLength, len - totalLength, newName);
    }

    private String correct(String rawCb) {
        if (whitelist.contains(rawCb))
            return rawCb;

        if (maxMismatch < 1)
            return null;

        char[] chars = rawCb.toCharArray();
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
