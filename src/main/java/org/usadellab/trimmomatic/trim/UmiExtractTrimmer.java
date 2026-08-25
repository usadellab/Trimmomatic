package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * UMIEXTRACT:<length>[:<separator>]
 *
 * Extracts a Unique Molecular Identifier (UMI) of <length> bases from the
 * 5' end of the read, appends it to the read name (separated by <separator>,
 * default "_"), and trims those bases from the read.
 *
 * The UMI is appended as:  @original_name<separator>UMI:<bases>
 *
 * For paired-end single-cell protocols (10x Genomics, Drop-seq, etc.) this
 * step should typically be applied only to R1, which contains the cell
 * barcode + UMI.  Apply it as one of the trimmer steps in the normal
 * pipeline; the read name change propagates to both R1 and R2 via downstream
 * tools that rely on matching names.
 *
 * Examples:
 *   UMIEXTRACT:12       – extract 12-base UMI, append as _UMI:ACGTACGTACGT
 *   UMIEXTRACT:10:__    – extract 10-base UMI, use __ as separator
 */
public class UmiExtractTrimmer extends AbstractSingleRecordTrimmer {
    private final int umiLength;
    private final String separator;

    public UmiExtractTrimmer(String args) {
        String[] arg = args.split(":");
        umiLength = Integer.parseInt(arg[0]);
        if (umiLength < 1)
            throw new IllegalArgumentException(
                    "UMIEXTRACT length must be >= 1, got: " + umiLength);
        separator = arg.length > 1 ? arg[1] : "_";
    }

    /**
     * UMIEXTRACT renames based on each record's own leading bases, so applying it
     * to both mates of a pair independently (as AbstractSingleRecordTrimmer would)
     * tags them with different, mismatched suffixes and desyncs mate names with
     * no error raised. It also chops real bases off whichever mate is not the
     * barcode/UMI read. This step therefore refuses. Route it to a single mate
     * via -pe1steps/-pe2steps or -technicalread, not via the shared list.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMIEXTRACT cannot run in symmetric paired-end mode: it would rename each mate "
                    + "independently from its own leading bases, desynchronising mate names. "
                    + "Route it to one mate only, e.g. via -pe1steps/-pe2steps or -technicalread.");
        return super.processRecords(in);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len <= umiLength)
            return null; // read too short to contain a UMI and any payload

        // Materialise only the UMI portion. getSequence() on the full view,
        // then substring just the UMI bases.
        String umi = in.getSequence().substring(0, umiLength);
        String newName = in.getName() + separator + "UMI:" + umi;

        // Return a view that skips the UMI bases, with the modified name.
        return new FastqRecord(in, umiLength, len - umiLength, newName);
    }
}
