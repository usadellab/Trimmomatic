package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Fringe-case tests for PolyXTrimmer.
 *
 * Covered here:
 *   - View records as input, including chained views
 *   - Run exactly minLength (trim) vs minLength-1 (no trim)
 *   - minLength=1 semantics
 *   - Single-character reads
 *   - Interior runs (no 3' tail) pass through unchanged
 *   - Run starting at position 1 (leaves a 1-base read)
 *   - Quality string is trimmed in sync with the sequence
 *   - Lowercase base in args is normalised to uppercase
 */
public class PolyXTrimmerEdgeCaseTest {

    private FastqRecord makeRecord(String seq, String qual) {
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    private FastqRecord makeRecord(String seq) {
        return makeRecord(seq, "I".repeat(seq.length()));
    }

    // ------------------------------------------------------------------
    // Boundary: exactly minLength vs minLength-1

    @Test
    public void testRunExactlyMinLengthTrims() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("GCGCAAAAA"); // exactly 5 A's
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    @Test
    public void testRunOneLessThanMinLengthPassesUnchanged() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("GCGCAAAA"); // 4 A's, one short
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGCAAAA", result.getSequence());
    }

    @Test
    public void testRunOneMoreThanMinLengthTrims() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("GCGCAAAAAA"); // 6 A's
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // minLength = 1

    @Test
    public void testMinLengthOneTrimsAnySingleTrailingBase() {
        PolyXTrimmer trimmer = new PolyXTrimmer("T:1");
        FastqRecord rec = makeRecord("ACGCT"); // 1 trailing T
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGC", result.getSequence());
    }

    @Test
    public void testMinLengthOneSingleTargetBaseDropped() {
        PolyXTrimmer trimmer = new PolyXTrimmer("T:1");
        FastqRecord rec = makeRecord("T"); // entire read is the target base
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testMinLengthOneSingleNonTargetPassesUnchanged() {
        PolyXTrimmer trimmer = new PolyXTrimmer("T:1");
        FastqRecord rec = makeRecord("A"); // not the target
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("A", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Interior-only runs (no 3' tail)

    @Test
    public void testInteriorRunIgnored() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("AAAAAACGT"); // run at 5' end, not 3'
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("AAAAAACGT", result.getSequence());
    }

    @Test
    public void testInteriorRunIgnoredMidRead() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("GCAAAAAGC"); // interior run, last char 'C'
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCAAAAAGC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Run starts at position 1 (leaves a 1-base prefix)

    @Test
    public void testRunStartsAtPositionOneLeavesSingleBase() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("CAAAAAA"); // 'C' then 6 A's
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result); // trimPos=1 ≠ 0 → not dropped
        assertEquals("C", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Quality string is trimmed in sync

    @Test
    public void testQualityTrimmedCorrectly() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:3");
        // "GCGCAAA", quality: "IIIIJJJ"
        // After trim: "GCGC" with quality "IIII"
        FastqRecord rec = makeRecord("GCGCAAA", "IIIIJJJ");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
        assertEquals("IIII", result.getQuality());
    }

    // ------------------------------------------------------------------
    // View records

    @Test
    public void testViewRecordTrimsTail() {
        // Base: 10 padding chars + "ACGTAAAAAA" (run of 6 A's)
        FastqRecord base = makeRecord("XXXXXXXXXXACGTAAAAAA"); // 20 chars
        FastqRecord view = new FastqRecord(base, 10, 10);      // "ACGTAAAAAA"
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord result = trimmer.processRecord(view);
        assertNotNull(result);
        assertEquals("ACGT", result.getSequence());
    }

    @Test
    public void testViewRecordNoTailPassesUnchanged() {
        FastqRecord base = makeRecord("XXXXXXXXXXACGTACGTAC");
        FastqRecord view = new FastqRecord(base, 10, 10); // "ACGTACGTAC"
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord result = trimmer.processRecord(view);
        assertNotNull(result);
        assertEquals("ACGTACGTAC", result.getSequence());
    }

    @Test
    public void testViewOfViewTrimsTail() {
        // Chained view: base → mid (skip 5) → leaf (skip 5 more) = chars 10-19
        FastqRecord base = makeRecord("XXXXXXXXXXXXXXXXXXTTTTTTTTTAAAA"); // 31 chars, 9 T's at end after offset
        // Actually: let me build it cleanly
        FastqRecord base2 = makeRecord("XXXXXX" + "GCGCGGGGGG"); // 16 chars total
        // View 1: skip 6, take 10 = "GCGCGGGGGG"
        FastqRecord mid  = new FastqRecord(base2, 6, 10);
        // View 2: skip 0, take 10 = same "GCGCGGGGGG"
        FastqRecord leaf = new FastqRecord(mid, 0, 10);
        PolyXTrimmer trimmer = new PolyXTrimmer("G:5");
        // "GCGCGGGGGG", 6 G's at end ≥ 5 → trim to "GCGC"
        FastqRecord result = trimmer.processRecord(leaf);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    @Test
    public void testViewEntireRunDropped() {
        // View contains only the target base
        FastqRecord base = makeRecord("XXXXXXXXXXXXXXXXAAAAAAAAAAA"); // 10 A's at end
        FastqRecord view = new FastqRecord(base, 17, 10); // "AAAAAAAAAA", wait, let me check chars
        // base = 17 X's + 10 A's = 27 chars  ...actually:
        FastqRecord base2 = makeRecord("XXXXXXXXXXXXXXXXX" + "AAAAAAAAAA"); // 17 X + 10 A = 27
        FastqRecord view2 = new FastqRecord(base2, 17, 10); // "AAAAAAAAAA"
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        assertNull(trimmer.processRecord(view2)); // entire view is poly-A → drop
    }

    // ------------------------------------------------------------------
    // PolyG (two-colour Illumina chemistry artefact)

    @Test
    public void testPolyGExactBoundaryTrimmed() {
        PolyXTrimmer trimmer = new PolyXTrimmer("G:10");
        FastqRecord rec = makeRecord("ACGTACGTACGGGGGGGGGG"); // exactly 10 G's
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGTAC", result.getSequence());
    }

    @Test
    public void testPolyGNineNotTrimmed() {
        PolyXTrimmer trimmer = new PolyXTrimmer("G:10");
        FastqRecord rec = makeRecord("ACGTACGTACGGGGGGGGG"); // 9 G's, one short
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGTACGGGGGGGGG", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Multiple bases at 3' end, target is interrupted

    @Test
    public void testRunInterruptedByOtherBase() {
        // "ACGTAAAGAA", last run of A's is 2 ('A','A' at end), not 5 total
        PolyXTrimmer trimmer = new PolyXTrimmer("A:3");
        FastqRecord rec = makeRecord("ACGTAAAGAA"); // trailing run = 2 A's < 3
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTAAAGAA", result.getSequence()); // unchanged
    }

    // ------------------------------------------------------------------
    // Lowercase base spec normalised

    @Test
    public void testLowercaseBaseArgNormalised() {
        PolyXTrimmer trimmer = new PolyXTrimmer("c:4"); // lowercase 'c'
        FastqRecord rec = makeRecord("ACGTCCCC"); // 4 trailing C's
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGT", result.getSequence());
    }
}
