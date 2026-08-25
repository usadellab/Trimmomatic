package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Quality-based trimmer using Trimmomatic's Maximum Information algorithm.
 *
 * The score arrays grow on demand when a read longer than the current capacity
 * is encountered, so reads of any length are handled safely.  Array growth is
 * thread-safe: the volatile {@code scoreArrays} field is read lock-free by
 * {@link #processRecord}, while {@link #grow} serialises rebuilds with a
 * {@code synchronized} lock and re-checks the size guard inside the lock to
 * prevent redundant rebuilds.
 */
public class MaximumInformationTrimmer extends AbstractSingleRecordTrimmer {

    // Initial length-score array capacity.  Grows automatically for longer reads.
    private static final int INITIAL_CAPACITY = 1000;

    public static final int MAXQUAL = 60;

    private final int parLength;
    private final float strictness;

    // Quality-probability weights: size MAXQUAL+1, fixed for the lifetime of this trimmer.
    private final double[] qualProbTmp;

    // Immutable snapshot of normalised long score arrays.  Published via a
    // volatile write so reader threads always see the latest grown version.
    private record ScoreArrays(long[] lengthScore, long[] qualProb) {}
    private volatile ScoreArrays scoreArrays;

    // ------------------------------------------------------------------
    // Normalisation helpers

    private static double calcNormalization(double[] array, int margin) {
        double maxVal = array[0];
        for (int i = 1; i < array.length; i++) {
            double val = Math.abs(array[i]);
            if (val > maxVal) {
                maxVal = val;
            }
        }
        return Long.MAX_VALUE / (maxVal * margin);
    }

    private static long[] normalize(double[] array, double ratio) {
        long[] out = new long[array.length];
        for (int i = 0; i < array.length; i++) {
            out[i] = (long) (array[i] * ratio);
        }
        return out;
    }

    // ------------------------------------------------------------------
    // Dynamic growth

    /**
     * Build a fresh pair of score arrays for the given read-length capacity.
     * The normalization margin is {@code capacity * 2} so that no long overflow
     * occurs when accumulating up to {@code capacity} scores.
     */
    private ScoreArrays buildArrays(int capacity) {
        double[] lengthScoreTmp = new double[capacity];
        for (int i = 0; i < capacity; i++) {
            // Unique weighting: logistic function on difference between length and parLength
            double pow1    = Math.exp(parLength - i - 1);
            double unique  = Math.log(1.0 / (1.0 + pow1));
            // Coverage weighting: length diluted by strictness
            double coverage = Math.log(i + 1) * (1 - strictness);
            lengthScoreTmp[i] = unique + coverage;
        }
        double normRatio = Math.max(
                calcNormalization(lengthScoreTmp, capacity * 2),
                calcNormalization(qualProbTmp,    capacity * 2));
        return new ScoreArrays(
                normalize(lengthScoreTmp, normRatio),
                normalize(qualProbTmp,    normRatio));
    }

    /**
     * Grow the score arrays to at least {@code needed} entries.
     * Synchronized so only one thread rebuilds at a time; the size guard is
     * re-checked inside the lock to avoid redundant rebuilds.
     */
    private synchronized void grow(int needed) {
        if (needed <= scoreArrays.lengthScore().length) {
            return; // another thread already grew it
        }
        int newCapacity = Math.max(needed, scoreArrays.lengthScore().length * 2);
        scoreArrays = buildArrays(newCapacity);
    }

    // ------------------------------------------------------------------
    // Constructor

    public MaximumInformationTrimmer(String args) {
        String[] arg = args.split(":");
        parLength  = Integer.parseInt(arg[0]);
        strictness = Float.parseFloat(arg[1]);

        qualProbTmp = new double[MAXQUAL + 1];
        for (int i = 0; i < qualProbTmp.length; i++) {
            // Quality weighting: probability of correctness, depending on strictness
            qualProbTmp[i] = Math.log(1 - Math.pow(0.1, (0.5 + i) / 10.0)) * strictness;
        }

        scoreArrays = buildArrays(INITIAL_CAPACITY);
    }

    // ------------------------------------------------------------------
    // Trimmer logic

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int[] quals = in.getQualityAsInteger(true);

        // Volatile read, no lock needed for the common (no-grow) path.
        ScoreArrays sa = scoreArrays;
        if (quals.length > sa.lengthScore().length) {
            grow(quals.length);
            sa = scoreArrays; // re-read after potential rebuild
        }

        long   accumQuality    = 0;
        double maxScore        = -Double.MAX_VALUE;
        int    maxScorePosition = 0;

        for (int i = 0; i < quals.length; i++) {
            int q = quals[i];
            if (q < 0) {
                q = 0;
            } else if (q > MAXQUAL) {
                q = MAXQUAL;
            }

            accumQuality += sa.qualProb()[q];
            long ls    = sa.lengthScore()[i];
            long score = ls + accumQuality;

            if (score >= maxScore) {
                maxScore        = score;
                maxScorePosition = i + 1;
            }
        }

        if (maxScorePosition < 1 || maxScore == 0.0) {
            return null;
        }
        if (maxScorePosition < quals.length) {
            return new FastqRecord(in, 0, maxScorePosition);
        }
        return in;
    }
}
