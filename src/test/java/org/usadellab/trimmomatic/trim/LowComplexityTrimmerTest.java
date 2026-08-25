package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class LowComplexityTrimmerTest {

    private FastqRecord makeRecord(String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    // ------------------------------------------------------------------
    // Homopolymer reads, entropy = 0.0

    @Test
    public void testHomopolymerDroppedWithThreshold1() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.0");
        FastqRecord rec = makeRecord("AAAAAAAAAAAAAAAA");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testHomopolymerDroppedWithThreshold0_5() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        FastqRecord rec = makeRecord("GGGGGGGGGGGGGGGG");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Di-nucleotide repeats, entropy = 1.0

    @Test
    public void testDiNucRepeatPassesAtExactThreshold() {
        // AT repeat has Shannon entropy exactly 1.0.
        // The comparison is >= so a read at the boundary passes.
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.0");
        FastqRecord rec = makeRecord("ATATATATATATATAT");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testDiNucRepeatDroppedAboveThreshold() {
        // AT repeat entropy 1.0 is below threshold 1.5 → dropped
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.5");
        FastqRecord rec = makeRecord("ATATATATATATATAT");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testDiNucRepeatPassesBelowThreshold() {
        // AT repeat: entropy 1.0 → passes when threshold = 0.5
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        FastqRecord rec = makeRecord("ATATATATATATATAT");
        assertNotNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Uniform ACGT, entropy = 2.0

    @Test
    public void testUniformAcgtAlwaysPasses() {
        // 16 bases, 4 of each, maximum entropy of 2.0
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.5");
        FastqRecord rec = makeRecord("ACGTACGTACGTACGT");
        assertNotNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // N-only reads

    @Test
    public void testAllNReadDropped() {
        // All Ns: total ACGT count = 0 → drop
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("NNNNNNNN");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Mixed N reads, Ns excluded from entropy calculation

    @Test
    public void testNsExcludedFromEntropy() {
        // "AAANNN": only A's among ACGT bases → entropy 0 → drop at threshold 0.5
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        FastqRecord rec = makeRecord("AAANNN");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Edge cases

    @Test
    public void testEmptyReadDropped() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testThresholdZeroPassesAnything() {
        // Even a homopolymer passes when minEntropy = 0.0
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("AAAAAAAAAA");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testThresholdTwoOnlyUniformPasses() {
        // Only a perfectly uniform ACGT read reaches entropy 2.0
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("2.0");
        FastqRecord rec = makeRecord("ACGT");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testThresholdTwoDropsNonUniform() {
        // 75% A, 25% C: entropy < 2.0 → drop
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("2.0");
        FastqRecord rec = makeRecord("AAAC");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Invalid args

    @Test
    public void testNegativeThresholdThrows() {
        assertThrows(IllegalArgumentException.class, () -> new LowComplexityTrimmer("-0.1"));
    }

    @Test
    public void testThresholdAboveTwoThrows() {
        assertThrows(IllegalArgumentException.class, () -> new LowComplexityTrimmer("2.1"));
    }
}
