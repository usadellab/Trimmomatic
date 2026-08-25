package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * UMISPLIT:<cbLength>:<umiLength>[:<separator>]
 *
 * Extracts a cell barcode of <cbLength> bases followed by a UMI of
 * <umiLength> bases, both from the 5' end of the read, appends both to the
 * read name (separated by <separator>, default "_"), and trims both off the
 * sequence. Unlike UMIEXTRACT, no leftover payload is required, a read
 * that is exactly cbLength + umiLength bases (the common case for droplet
 * platforms where the whole read IS the barcode/UMI construct, e.g. 10x
 * Chromium R1) survives with a zero-length sequence; it is dropped only if
 * shorter than that.
 *
 * This is a separate step from UMIEXTRACT. An overloaded 2-arg form was not
 * possible, because UMIEXTRACT:<length>:<separator> already uses its 2nd
 * colon-arg for a separator string. Reusing that shape for a second length
 * would be ambiguous to the parser and to a human reading a script.
 *
 * Does not whitelist-correct the barcode, see UMIDROPLETCORRECT for that.
 *
 * The read is renamed as: @original_name<separator><CB bases><separator><UMI bases>
 * These are bare sequences, no "CB:"/"UMI:" labels, matching umi_tools extract's own
 * convention so umi_tools dedup/count's default read-name parser (last
 * underscore-delimited, expects a bare nucleotide string) can consume this
 * directly.
 *
 * Examples:
 *   UMISPLIT:16:12       - 16bp CB + 12bp UMI, appended as _<16bp>_<12bp>
 *   UMISPLIT:12:8:__     - 12bp CB + 8bp UMI (Drop-seq style), __ as separator
 */
public class UmiSplitTrimmer extends AbstractSingleRecordTrimmer {
    private final int cbLength;
    private final int umiLength;
    private final String separator;

    public UmiSplitTrimmer(String args) {
        String[] arg = args.split(":");
        if (arg.length < 2)
            throw new IllegalArgumentException(
                    "UMISPLIT requires two lengths: <cbLength>:<umiLength>, got: " + args);

        cbLength = Integer.parseInt(arg[0]);
        umiLength = Integer.parseInt(arg[1]);

        if (cbLength < 1)
            throw new IllegalArgumentException("UMISPLIT cbLength must be >= 1, got: " + cbLength);
        if (umiLength < 1)
            throw new IllegalArgumentException("UMISPLIT umiLength must be >= 1, got: " + umiLength);

        separator = arg.length > 2 ? arg[2] : "_";
    }

    /**
     * Same symmetric-mode hazard as UMIEXTRACT, see its processRecords() for
     * the full explanation. Refusing keeps mate names in sync.
     */
    @Override
    public FastqRecord[] processRecords(FastqRecord[] in) {
        if (in != null && in.length > 1)
            throw new IllegalStateException(
                    "UMISPLIT cannot run in symmetric paired-end mode: it would rename each mate "
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
        String cb = seq.substring(0, cbLength);
        String umi = seq.substring(cbLength, totalLength);
        // Bare sequences, no "CB:"/"UMI:" labels, matches umi_tools extract's own
        // convention (@ReadName_BARCODE_UMI) so umi_tools dedup/count's default
        // read-name parser (last-underscore-delimited, expects a bare nucleotide
        // string) can consume this directly without a label getting in the way.
        String newName = in.getName() + separator + cb + separator + umi;

        return new FastqRecord(in, totalLength, len - totalLength, newName);
    }
}
