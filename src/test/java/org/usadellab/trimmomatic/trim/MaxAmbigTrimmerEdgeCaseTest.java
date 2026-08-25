package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Fringe-case tests for MaxAmbigTrimmer.
 *
 * Covered here:
 *   - View records (trimmed FastqRecord), ensures getLength() not sequence.length()
 *   - Floating-point precision at the boundary (1/3 vs 0.333 / 0.334)
 *   - Exact fraction equality uses <= (passes at exact boundary)
 *   - Case sensitivity: lowercase 'n' is NOT counted, only uppercase 'N'
 *   - Chain-of-views (view of a view)
 *   - Single-base reads
 *   - Very large reads
 *   - All-non-ACGTN characters (no N counted)
 */
public class MaxAmbigTrimmerEdgeCaseTest {

    private FastqRecord makeRecord(String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    // ------------------------------------------------------------------
    // View record behaviour

    @Test
    public void testViewRecordNoNInViewButNsInBase() {
        // Base has 10 N's at the end; the view covers only the clean prefix.
        FastqRecord base = makeRecord("ACGTACGTACNNNNNNNNNN"); // 20 chars, 10 N's
        FastqRecord view = new FastqRecord(base, 0, 10);      // "ACGTACGTAC", 0 N's
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        // View has 0 N's → fraction 0.0 ≤ 0.0 → pass
        assertNotNull(trimmer.processRecord(view));
    }

    @Test
    public void testViewRecordAllNInView() {
        // Base: first 10 chars clean, last 10 all N.
        FastqRecord base = makeRecord("ACGTACGTACNNNNNNNNNN");
        FastqRecord view = new FastqRecord(base, 10, 10); // "NNNNNNNNNN", 10/10 N's
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.5");
        // 10 N's / 10 bases = 1.0 > 0.5 → drop
        assertNull(trimmer.processRecord(view));
    }

    @Test
    public void testViewRecordFractionWithinThreshold() {
        // Base: "ACGTNN" + 14 clean bases → 20 chars
        FastqRecord base = makeRecord("ACGTNNAAAAAAAAAAAAAAA"); // 21 chars, 2 Ns at pos 4-5
        // View: chars 0-9 → "ACGTNNAAAA", 2 N's out of 10 = 0.2
        FastqRecord view = new FastqRecord(base, 0, 10);
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        // 2/10 = 0.2 ≤ 0.2 → pass (boundary inclusive)
        assertNotNull(trimmer.processRecord(view));
    }

    @Test
    public void testChainedViewBothClean() {
        // View of a view, ensure no offset arithmetic bugs.
        FastqRecord base = makeRecord("NNNNNNNNNNACGTACGTAC"); // 20 chars
        FastqRecord mid  = new FastqRecord(base, 5, 15);       // "NNNNNACGTACGTAC", 5/15 N's
        FastqRecord leaf = new FastqRecord(mid, 5, 10);        // "ACGTACGTAC", 0/10 N's
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        assertNotNull(trimmer.processRecord(leaf));
    }

    // ------------------------------------------------------------------
    // Floating-point boundary precision

    @Test
    public void testOneNInThreeBases_belowThresholdDrops() {
        // 1/3 = 0.333...  The trimmer uses <=.
        // 0.333... > 0.333 → drop.
        FastqRecord rec = makeRecord("ACN");
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.333");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testOneNInThreeBases_aboveThresholdPasses() {
        // 1/3 = 0.333... ≤ 0.334 → pass.
        FastqRecord rec = makeRecord("ACN");
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.334");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testExactFractionEqualityPassesBoundary() {
        // 2 N's / 10 bases = exactly 0.2.  <= comparison → pass.
        FastqRecord rec = makeRecord("ACGTACGNNN"); // 3 Ns actually... let me use 2 Ns in 10
        // "ACGTACGTNNAAAAAA" wait, let me keep it simple: 2 N's in 10 chars
        FastqRecord rec2 = makeRecord("ACGTACGNN" + "A"); // "ACGTACGNNA" = 2 Ns / 10
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        assertNotNull(trimmer.processRecord(rec2));
    }

    @Test
    public void testOneAboveBoundaryDrops() {
        // 3 N's out of 10 = 0.3 > 0.2 → drop
        FastqRecord rec = makeRecord("ACGTNNNACG"); // 3 N's / 10
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Case sensitivity

    @Test
    public void testLowercaseNNotCounted() {
        // Lowercase 'n' should NOT be treated as ambiguous by the trimmer.
        // Only uppercase 'N' is counted.
        FastqRecord rec = makeRecord("ACGTACGTnnnnn"); // 5 lowercase n's, 0 uppercase N's
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        // 0 uppercase N's / 13 bases = 0.0 ≤ 0.0 → pass
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testMixedCaseNOnlyUpperCounted() {
        // 2 uppercase N's and 3 lowercase n's in a 10-base read
        FastqRecord rec = makeRecord("ACNTnNGCAT"); // N at 2, n at 4, N at 5 → 2 uppercase N's
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        // 2/10 = 0.2 ≤ 0.2 → pass
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testMixedCaseStrictlyUpperCounted_drops() {
        // 3 uppercase N's / 10 bases = 0.3 > 0.2 → drop
        FastqRecord rec = makeRecord("NNNACGTnnn"); // 3 uppercase Ns
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Single-base reads

    @Test
    public void testSingleN_zeroThreshold_drops() {
        FastqRecord rec = makeRecord("N"); // 1/1 = 1.0 > 0.0
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testSingleA_zeroThreshold_passes() {
        FastqRecord rec = makeRecord("A"); // 0/1 = 0.0 ≤ 0.0
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testSingleN_fullThreshold_passes() {
        FastqRecord rec = makeRecord("N"); // 1/1 = 1.0 ≤ 1.0
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("1.0");
        assertNotNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Large read, very low N fraction

    @Test
    public void testOneNInHundredBases() {
        String seq = "A".repeat(99) + "N"; // 1/100 = 0.01
        FastqRecord rec = makeRecord(seq);
        MaxAmbigTrimmer trimmer009 = new MaxAmbigTrimmer("0.009");
        MaxAmbigTrimmer trimmer010 = new MaxAmbigTrimmer("0.01");
        assertNull(trimmer009.processRecord(rec));     // 0.01 > 0.009 → drop
        assertNotNull(trimmer010.processRecord(rec));  // 0.01 ≤ 0.01  → pass
    }

    // ------------------------------------------------------------------
    // Non-ACGTN characters are never counted

    @Test
    public void testNonStandardBases_notCountedAsN() {
        // IUPAC ambiguity codes like 'R','Y','S','W' are not 'N' → not counted
        FastqRecord rec = makeRecord("ACGTRYSWKM"); // no uppercase N → fraction 0.0
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        assertNotNull(trimmer.processRecord(rec));
    }
}
