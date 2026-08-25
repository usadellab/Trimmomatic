package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Fringe-case tests for LowComplexityTrimmer.
 *
 * Covered here:
 *   - View records as input
 *   - Shannon entropy boundary arithmetic for 2-, 3-, and 4-base alphabets
 *   - Entropy is inclusive at the threshold (entropy >= threshold passes)
 *   - Reads with only 2 unequal base proportions
 *   - Reads with 3 equal bases (entropy = log2(3) ≈ 1.585)
 *   - Reads with 4 equal bases (entropy = 2.0)
 *   - Single-base reads (entropy = 0.0)
 *   - Reads where all counted bases are the same after N exclusion
 *   - Chained view whose low/high complexity parts are tested independently
 */
public class LowComplexityTrimmerEdgeCaseTest {

    private FastqRecord makeRecord(String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    // ------------------------------------------------------------------
    // Single-base reads (entropy = 0.0)

    @Test
    public void testSingleBaseRead_entropy0_thresholdZero_passes() {
        // entropy 0.0 >= 0.0 → pass
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("A");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testSingleBaseRead_entropy0_positiveThreshold_drops() {
        // entropy 0.0 < 0.5 → drop
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        FastqRecord rec = makeRecord("C");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // 2-base alphabet

    @Test
    public void testTwoEqualBases_entropy1_exactThreshold_passes() {
        // 8 A's + 8 C's = entropy exactly 1.0; threshold 1.0 → 1.0 >= 1.0 → pass
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.0");
        FastqRecord rec = makeRecord("AAAAAAAACCCCCCCC");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testTwoEqualBases_entropy1_aboveThreshold_drops() {
        // entropy 1.0 < 1.1 → drop
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.1");
        FastqRecord rec = makeRecord("AAAAAAAACCCCCCCC");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testTwoUnequalBases_computedEntropy() {
        // 12 A's + 4 C's (total=16): p_A=0.75, p_C=0.25
        // H = -(0.75*log2(0.75) + 0.25*log2(0.25))
        //   = -(0.75*(-0.415) + 0.25*(-2))
        //   = -(-0.311 + -0.5) = 0.811
        LowComplexityTrimmer trimmerLow  = new LowComplexityTrimmer("0.8");
        LowComplexityTrimmer trimmerHigh = new LowComplexityTrimmer("0.82");
        FastqRecord rec = makeRecord("AAAAAAAAAAAACCCC"); // 12 A, 4 C
        assertNotNull(trimmerLow.processRecord(rec));  // 0.811 >= 0.80 → pass
        assertNull(trimmerHigh.processRecord(rec));   // 0.811 <  0.82 → drop
    }

    // ------------------------------------------------------------------
    // 3-base alphabet, entropy = log2(3) ≈ 1.5849625

    @Test
    public void testThreeEqualBases_belowLog2of3_passes() {
        // 4 A's + 4 C's + 4 G's (T=0), entropy = log2(3) ≈ 1.585
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.58");
        FastqRecord rec = makeRecord("AAAACCCCGGGG");
        assertNotNull(trimmer.processRecord(rec));  // 1.585 >= 1.58 → pass
    }

    @Test
    public void testThreeEqualBases_aboveLog2of3_drops() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.59");
        FastqRecord rec = makeRecord("AAAACCCCGGGG");
        assertNull(trimmer.processRecord(rec));     // 1.585 < 1.59 → drop
    }

    // ------------------------------------------------------------------
    // 4-base alphabet, entropy = 2.0 exactly

    @Test
    public void testFourEqualBases_entropy2_passes() {
        // 4 of each: A, C, G, T → entropy = 2.0
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("2.0");
        FastqRecord rec = makeRecord("AAAACCCCGGGGTTTT");
        assertNotNull(trimmer.processRecord(rec));  // 2.0 >= 2.0 → pass
    }

    @Test
    public void testFourUnequalBases_entropy1_75() {
        // 8A, 4C, 2G, 2T (total=16)
        // p_A=0.5, p_C=0.25, p_G=0.125, p_T=0.125
        // H = 0.5+0.5+0.375+0.375 = 1.75
        LowComplexityTrimmer trimmerBelow = new LowComplexityTrimmer("1.74");
        LowComplexityTrimmer trimmerAbove = new LowComplexityTrimmer("1.76");
        FastqRecord rec = makeRecord("AAAAAAAACCCCGGTT"); // 8A 4C 2G 2T
        assertNotNull(trimmerBelow.processRecord(rec)); // 1.75 >= 1.74 → pass
        assertNull(trimmerAbove.processRecord(rec));   // 1.75 <  1.76 → drop
    }

    // ------------------------------------------------------------------
    // Reads with N's, N bases excluded from entropy

    @Test
    public void testReadWithNs_onlyClearedBasesCountedForEntropy() {
        // "AAAANNNN" → countA=4, total=4 → entropy=0.0 → drop at threshold 0.5
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        FastqRecord rec = makeRecord("AAAANNNN");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testReadMostlyNs_highEntropyBases_passes() {
        // 4 each of A, C, G, T buried in N's → entropy 2.0 from the 16 non-N bases
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.9");
        FastqRecord rec = makeRecord("NNNNNNNNNNNACGTACGTACGTACGT"); // many N's
        assertNotNull(trimmer.processRecord(rec)); // entropy from ACGT = 2.0 → pass
    }

    // ------------------------------------------------------------------
    // View records

    @Test
    public void testViewRecord_lowComplexityPart_dropped() {
        // Base: 10 clean chars + 10 all-A chars
        FastqRecord base = makeRecord("ACGTACGTACAAAAAAAAAA"); // 20 chars
        FastqRecord viewLow = new FastqRecord(base, 10, 10); // "AAAAAAAAAA", entropy 0
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.5");
        assertNull(trimmer.processRecord(viewLow));
    }

    @Test
    public void testViewRecord_highComplexityPart_passes() {
        // Base: 10 all-A + 10 uniform ACGT
        FastqRecord base = makeRecord("AAAAAAAAAACGTACGTACG"); // 20 chars
        FastqRecord viewHigh = new FastqRecord(base, 10, 10); // "CGTACGTACG", entropy ~2
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.5");
        assertNotNull(trimmer.processRecord(viewHigh));
    }

    @Test
    public void testViewOfView_correctEntropyUsed() {
        // base: 5 junk + 4A + 4C + 4G + 3T (= 15 non-junk chars)
        // We take the middle 12 chars (all 4A + 4C + 4G = entropy log2(3) ≈ 1.585)
        FastqRecord base = makeRecord("XXXXXAAAACCCCGGGGTT" + "T"); // 21 chars
        FastqRecord mid  = new FastqRecord(base, 5, 12); // "AAAACCCCGGGG"
        FastqRecord leaf = new FastqRecord(mid, 0, 12);  // same
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.58");
        assertNotNull(trimmer.processRecord(leaf)); // 1.585 >= 1.58 → pass
    }

    // ------------------------------------------------------------------
    // Reads with only non-ACGT bases (all N or all ambiguity codes)

    @Test
    public void testAllNRead_dropsAtAnyPositiveThreshold() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("NNNNNNNNN");
        // total ACGT = 0 → the all-N special-case returns null
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testOnlyAmbiguityCodes_noACGT_dropped() {
        // No standard bases at all
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("0.0");
        FastqRecord rec = makeRecord("RYSWKM"); // IUPAC ambiguity, no A/C/G/T
        // countA=countC=countG=countT=0 → total=0 → drop
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Two-base equal alphabet via GC-only read

    @Test
    public void testGCOnlyRead_entropy1_exactThreshold_passes() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.0");
        FastqRecord rec = makeRecord("GCGCGCGCGCGCGCGC"); // equal G and C → entropy 1.0
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testATOnlyRead_entropy1_aboveThreshold_drops() {
        LowComplexityTrimmer trimmer = new LowComplexityTrimmer("1.5");
        FastqRecord rec = makeRecord("ATATATATATAT"); // equal A and T → entropy 1.0 < 1.5
        assertNull(trimmer.processRecord(rec));
    }
}
