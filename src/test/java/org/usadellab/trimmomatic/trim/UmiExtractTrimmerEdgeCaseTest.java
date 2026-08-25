package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Fringe-case tests for UmiExtractTrimmer.
 *
 * Covered here:
 *   - View records as input (including chained views), verifies view-offset arithmetic
 *   - UMI length = 1 (minimum useful value)
 *   - UMI length = read length (equals → drop)
 *   - UMI length = read length - 1 (leaves exactly 1-base payload)
 *   - Empty separator via trailing colon in args ("4:"), Java split semantics
 *   - Read name that already contains the separator string
 *   - UMI bases are all N (ambiguous)
 *   - headPos is accumulated correctly after extraction from a view
 *   - Quality string trimmed to the payload portion
 *   - Very long separator string
 */
public class UmiExtractTrimmerEdgeCaseTest {

    private FastqRecord makeRecord(String name, String seq, String qual) {
        return new FastqRecord(name, seq, "", qual, 33);
    }

    private FastqRecord makeRecord(String name, String seq) {
        return makeRecord(name, seq, "I".repeat(seq.length()));
    }

    // ------------------------------------------------------------------
    // UMI length edge values

    @Test
    public void testUmiLengthOne_extractsSingleBase() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("1");
        FastqRecord rec = makeRecord("r1", "ACGTACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("r1_UMI:A", result.getName());
        assertEquals("CGTACGT", result.getSequence());
        assertEquals(7, result.getLength());
    }

    @Test
    public void testUmiLengthEqualsReadLength_drops() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("8");
        FastqRecord rec = makeRecord("r1", "ACGTACGT"); // length = umiLength
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testUmiLengthReadLengthMinusOne_leavesOneBase() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("7");
        FastqRecord rec = makeRecord("r1", "ACGTACGT"); // 8 bases; UMI=7, payload=1
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("r1_UMI:ACGTACG", result.getName());
        assertEquals("T", result.getSequence());
        assertEquals(1, result.getLength());
    }

    @Test
    public void testUmiLongerThanRead_drops() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("12");
        FastqRecord rec = makeRecord("r1", "ACGTACGT"); // 8 < 12
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Separator edge cases

    @Test
    public void testDefaultSeparatorUnderscore() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord rec = makeRecord("read1", "TTTGACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertTrue(result.getName().startsWith("read1_UMI:"), result.getName());
    }

    @Test
    public void testCustomSeparatorDoubleUnderscore() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4:__");
        FastqRecord rec = makeRecord("read1", "TTTGACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("read1__UMI:TTTG", result.getName());
    }

    @Test
    public void testTrailingColonInArgs_defaultSeparatorUsed() {
        // "4:" → Java split(":") drops trailing empty tokens → arg.length == 1
        // → separator defaults to "_"
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4:");
        FastqRecord rec = makeRecord("r1", "ACGTGGGG");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        // separator should be "_" (default) because "4:".split(":") = ["4"]
        assertEquals("r1_UMI:ACGT", result.getName());
    }

    @Test
    public void testReadNameAlreadyContainsSeparator() {
        // Name with "_" already; separator is "_"; result should append, not replace
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord rec = makeRecord("read_1_extra", "TTTTGGGG");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("read_1_extra_UMI:TTTT", result.getName());
    }

    @Test
    public void testVeryLongSeparator() {
        // Separator must not contain colons, the constructor uses split(":") so
        // colons inside the separator would be treated as extra split points.
        String sep = "___LONG_SEP___";
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("3:" + sep);
        FastqRecord rec = makeRecord("r1", "NNNACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("r1" + sep + "UMI:NNN", result.getName());
        assertEquals("ACGT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // UMI content

    @Test
    public void testUmiAllN_appendedCorrectly() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord rec = makeRecord("r1", "NNNNACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("r1_UMI:NNNN", result.getName());
        assertEquals("ACGT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Quality string is trimmed to payload portion

    @Test
    public void testQualityTrimmedToPayload() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        String seq  = "UUUUACGTACGT"; // 12 bases (U = placeholder)
        String qual = "IIII" + "JJJJJJJJ"; // first 4 = UMI quality
        FastqRecord rec = makeRecord("r1", seq, qual);
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("JJJJJJJJ", result.getQuality());
    }

    // ------------------------------------------------------------------
    // View records, offset arithmetic

    @Test
    public void testViewRecord_umiExtractedFromViewStart() {
        // Base:  "AAAAAGGGGGCCCCCTTTTTAAAAAGGGGG" (30 chars)
        // View:   chars 10-29 = "CCCCCTTTTTAAAAAGGGGG"
        // UMI=5: extracts "CCCCC", payload = "TTTTTAAAAAGGGGG"
        String baseSeq = "AAAAAGGGGGCCCCCTTTTTAAAAAGGGGG";
        FastqRecord base = makeRecord("r1", baseSeq, "I".repeat(30));
        FastqRecord view = new FastqRecord(base, 10, 20);

        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("5");
        FastqRecord result = trimmer.processRecord(view);

        assertNotNull(result);
        assertEquals("r1_UMI:CCCCC", result.getName());
        assertEquals("TTTTTAAAAAGGGGG", result.getSequence());
        assertEquals(15, result.getLength());
    }

    @Test
    public void testChainedView_umiExtractedFromLeafStart() {
        // Chained views: base → skip 5 → skip 5 more → leaf starts at offset 10 of base
        // Same expected result as above
        String baseSeq = "AAAAAGGGGGCCCCCTTTTTAAAAAGGGGG";
        FastqRecord base = makeRecord("r1", baseSeq, "I".repeat(30));
        FastqRecord mid  = new FastqRecord(base, 5, 25);  // chars 5-29 = "GGGGGCCCCCTTTTTAAAAAGGGGG"
        FastqRecord leaf = new FastqRecord(mid, 5, 20);   // chars 10-29 = "CCCCCTTTTTAAAAAGGGGG"

        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("5");
        FastqRecord result = trimmer.processRecord(leaf);

        assertNotNull(result);
        assertEquals("r1_UMI:CCCCC", result.getName());
        assertEquals("TTTTTAAAAAGGGGG", result.getSequence());
    }

    @Test
    public void testViewRecord_qualityAlsoOffset() {
        // Verify the quality string of the result starts at the correct raw position.
        // Base quality: 30 chars of alternating 'I' and 'J'
        StringBuilder qb = new StringBuilder(30);
        for (int i = 0; i < 30; i++) qb.append(i % 2 == 0 ? 'I' : 'J');
        String baseQual = qb.toString(); // "IJIJIJ..."

        FastqRecord base = makeRecord("r1", "AAAAAGGGGGCCCCCTTTTTAAAAAGGGGG", baseQual);
        // View of chars 10-29 → rawQual substring: chars 10-29 of baseQual
        FastqRecord view = new FastqRecord(base, 10, 20);

        // view quality = baseQual.substring(10, 30)
        String expectedViewQual = baseQual.substring(10, 30);
        // After UMI extraction of 4 chars: payload quality = chars 4-19 of viewQual
        String expectedPayloadQual = expectedViewQual.substring(4);

        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord result = trimmer.processRecord(view);
        assertNotNull(result);
        assertEquals(expectedPayloadQual, result.getQuality());
    }

    @Test
    public void testViewRecord_dropsWhenPayloadEmpty() {
        // View of exactly umiLength bases → no payload → drop
        FastqRecord base = makeRecord("r1", "XXXXXXXXXXXXXXXXXACGTACGT"); // 25 chars
        FastqRecord view = new FastqRecord(base, 17, 8); // "ACGTACGT", 8 chars
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("8"); // UMI = 8 = view length
        assertNull(trimmer.processRecord(view));
    }

    // ------------------------------------------------------------------
    // headPos accumulation

    @Test
    public void testHeadPosAccumulatesAcrossUmiExtraction() {
        // The result's headPos should be base.headPos + umiLength.
        FastqRecord base = makeRecord("r1", "ACGTACGT");
        // First create a trimmed view with headPos=2 (trimmed 2 from start)
        FastqRecord view = new FastqRecord(base, 2, 6); // headPos=2, seq="GTACGT"
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("3"); // extract 3 from view
        FastqRecord result = trimmer.processRecord(view);
        assertNotNull(result);
        // headPos should be 2 (view offset) + 3 (UMI) = 5
        assertEquals(5, result.getHeadPos());
        assertEquals("CGT", result.getSequence()); // "GTACGT" without first 3 = "CGT"
    }
}
