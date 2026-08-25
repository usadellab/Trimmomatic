package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * POLYX:<base>:<minLength>
 *
 * Trims a homopolymer run of <base> from the 3' end of the read when the run
 * is at least <minLength> bases long.  If the entire read consists of that
 * base the read is dropped.
 *
 * Examples:
 *   POLYX:A:10  – trim polyA tails of ≥10 bases (common in RNA-seq / scRNA-seq)
 *   POLYX:T:8   – trim polyT artefacts from reverse-transcription primers
 *   POLYX:G:10  – trim polyG artefacts from two-colour Illumina chemistry
 */
public class PolyXTrimmer extends AbstractSingleRecordTrimmer {
    private final char base;
    private final int minLength;

    public PolyXTrimmer(String args) {
        String[] arg = args.split(":");
        if (arg.length < 2)
            throw new IllegalArgumentException(
                    "POLYX requires two arguments: <base>:<minLength>, got: " + args);

        String baseStr = arg[0].toUpperCase();
        if (baseStr.length() != 1 || "ACGTN".indexOf(baseStr.charAt(0)) < 0)
            throw new IllegalArgumentException(
                    "POLYX base must be a single nucleotide character (A/C/G/T/N), got: " + arg[0]);

        base = baseStr.charAt(0);
        minLength = Integer.parseInt(arg[1]);
        if (minLength < 1)
            throw new IllegalArgumentException(
                    "POLYX minLength must be >= 1, got: " + minLength);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len == 0)
            return null;

        String seq = in.getSequence();

        // Count the run of <base> at the 3' end
        int trimPos = len;
        while (trimPos > 0 && seq.charAt(trimPos - 1) == base)
            trimPos--;

        int runLength = len - trimPos;
        if (runLength < minLength)
            return in; // run too short, pass through unchanged

        if (trimPos == 0)
            return null; // entire read is poly-X, drop it

        return new FastqRecord(in, 0, trimPos);
    }
}
