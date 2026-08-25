package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * LOWCOMPLEXITY:<minEntropy>
 *
 * Drops the read if its Shannon entropy (calculated over A/C/G/T base
 * frequencies) falls below minEntropy.  N bases are excluded from the
 * calculation.
 *
 * Entropy ranges from 0.0 (single-base homopolymer, e.g. AAAAA…) to
 * 2.0 (perfectly uniform A/C/G/T distribution).
 *
 * Recommended thresholds:
 *   minEntropy=1.5  – aggressive; drops homopolymers AND simple di-nucleotide
 *                      repeats (ATATATATAT has entropy exactly 1.0)
 *   minEntropy=1.0  – moderate; drops homopolymers but passes di-nucleotide
 *                      repeats (boundary is inclusive: entropy >= minEntropy passes)
 *   minEntropy=0.5  – lenient; catches only near-total-homopolymer reads
 *
 * Example: LOWCOMPLEXITY:1.0
 */
public class LowComplexityTrimmer extends AbstractSingleRecordTrimmer {
    private final double minEntropy;

    public LowComplexityTrimmer(String args) {
        minEntropy = Double.parseDouble(args);
        if (minEntropy < 0.0 || minEntropy > 2.0)
            throw new IllegalArgumentException(
                    "LOWCOMPLEXITY minEntropy must be between 0.0 and 2.0, got: " + minEntropy);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len == 0)
            return null;

        String seq = in.getSequence();

        int countA = 0, countC = 0, countG = 0, countT = 0;
        for (int i = 0; i < len; i++) {
            switch (seq.charAt(i)) {
                case 'A' -> countA++;
                case 'C' -> countC++;
                case 'G' -> countG++;
                case 'T' -> countT++;
                // N and other ambiguity codes are excluded from entropy calculation
            }
        }

        int total = countA + countC + countG + countT;
        if (total == 0)
            return null; // read is all-N, drop it

        double entropy = shannonEntropy(total, countA, countC, countG, countT);

        return entropy >= minEntropy ? in : null;
    }

    private static double shannonEntropy(int total, int... counts) {
        double entropy = 0.0;
        for (int c : counts) {
            if (c > 0) {
                double p = (double) c / total;
                entropy -= p * (Math.log(p) / Math.log(2));
            }
        }
        return entropy;
    }
}
