package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Long-read safety tests for MaximumInformationTrimmer.
 *
 * Prior to V0.42 the trimmer pre-allocated a fixed array of 1000 entries and
 * crashed with ArrayIndexOutOfBoundsException for any read longer than 1000 bp.
 * These tests verify the dynamic-growth fix introduced in V0.42.
 *
 * Covered here:
 *   - Read exactly at initial capacity (1000 bp), no crash
 *   - Read one beyond initial capacity (1001 bp), triggers grow, no crash
 *   - Very long reads (5000 bp, 10 000 bp, 100 000 bp), no crash
 *   - Sequential reads of increasing length, all passes
 *   - Short read processed after a grow, still works
 *   - Concurrent access from multiple threads, no exception
 */
public class MaximumInformationTrimmerLongReadTest {

    /** Build a FastqRecord of exactly {@code length} bases at uniform Phred quality {@code q}. */
    private static FastqRecord makeRecord(int length, int q) {
        String seq  = "A".repeat(length);
        String qual = String.valueOf((char) (33 + q)).repeat(length);
        return new FastqRecord("r1", seq, "", qual, 33);
    }

    // ------------------------------------------------------------------
    // At and around initial capacity (1000 bp)

    @Test
    public void testReadExactlyInitialCapacity_noCrash() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        FastqRecord result = trimmer.processRecord(makeRecord(1000, 30));
        assertNotNull(result);
    }

    @Test
    public void testReadOneOverInitialCapacity_growsAndProcesses() {
        // Prior to the fix this threw ArrayIndexOutOfBoundsException.
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        FastqRecord result = trimmer.processRecord(makeRecord(1001, 30));
        assertNotNull(result);
    }

    @Test
    public void testReadTwoOverInitialCapacity_noCrash() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        assertNotNull(trimmer.processRecord(makeRecord(1002, 30)));
    }

    // ------------------------------------------------------------------
    // Very long reads

    @Test
    public void testVeryLongRead_5000bp_noCrash() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        assertNotNull(trimmer.processRecord(makeRecord(5000, 30)));
    }

    @Test
    public void testVeryLongRead_10000bp_noCrash() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("500:0.5");
        assertNotNull(trimmer.processRecord(makeRecord(10_000, 30)));
    }

    @Test
    public void testVeryLongRead_100kbp_noCrash() {
        // Simulates a full-length ONT read (~100 kb).
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("1000:0.5");
        // Must not crash regardless of trim result (may return null for very low-quality).
        trimmer.processRecord(makeRecord(100_000, 15));
    }

    // ------------------------------------------------------------------
    // Sequential grows

    @Test
    public void testSequentialGrowth_increasingLengths() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        for (int len : new int[]{ 800, 1100, 2200, 4500, 9100 }) {
            FastqRecord result = trimmer.processRecord(makeRecord(len, 30));
            assertNotNull(result, "Unexpected null for read length " + len);
        }
    }

    // ------------------------------------------------------------------
    // Short reads still work after a grow

    @Test
    public void testShortReadAfterGrowth_stillWorks() {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        trimmer.processRecord(makeRecord(2000, 30)); // force a grow
        assertNotNull(trimmer.processRecord(makeRecord(50, 30)));
    }

    // ------------------------------------------------------------------
    // Thread-safety: concurrent access to a single trimmer instance

    @Test
    public void testConcurrentAccess_noException() throws InterruptedException {
        MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("36:0.5");
        int nThreads = 8;
        Thread[] threads = new Thread[nThreads];
        Throwable[] errors = new Throwable[nThreads];

        for (int t = 0; t < nThreads; t++) {
            final int tid = t;
            threads[t] = Thread.ofVirtual().start(() -> {
                try {
                    // Mix of lengths that straddle the initial capacity to provoke concurrent grows.
                    for (int len : new int[]{ 900, 1050, 2000, 500, 3000 }) {
                        trimmer.processRecord(makeRecord(len, 30));
                    }
                } catch (Throwable e) {
                    errors[tid] = e;
                }
            });
        }

        for (Thread th : threads) {
            th.join(5000);
        }

        for (int t = 0; t < nThreads; t++) {
            if (errors[t] != null) {
                throw new AssertionError("Thread " + t + " threw: " + errors[t], errors[t]);
            }
        }
    }
}
