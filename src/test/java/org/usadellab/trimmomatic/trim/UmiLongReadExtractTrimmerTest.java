package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class UmiLongReadExtractTrimmerTest {

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    // ------------------------------------------------------------------
    // NONE anchor (undisclosed cassette, e.g. PCS114/PCB114)

    @Test
    public void testNoAnchor_extractsFromPositionZero() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNNNNNN:0:0");

        FastqRecord rec = makeRecord("r1", "ACGTACGT" + "PAYLOADSEQUENCE");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:ACGTACGT", result.getName());
        assertEquals("PAYLOADSEQUENCE", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Anchor search, exact position

    @Test
    public void testAnchorExact_extractsUmiAfterAnchor() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("GTATCGTGT:NNNN:0:0");

        FastqRecord rec = makeRecord("r1", "GTATCGTGT" + "AACC" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:AACC", result.getName());
        assertEquals("PAYLOAD", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Indel-shift cascade

    @Test
    public void testAnchorShiftedByOneBase_foundWithIndelShift() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("GTATCGTGT:NNNN:0:2");

        // One extra base inserted before the anchor -- anchor now starts at offset 1.
        FastqRecord rec = makeRecord("r1", "X" + "GTATCGTGT" + "AACC" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:AACC", result.getName());
        assertEquals("PAYLOAD", result.getSequence());
    }

    @Test
    public void testAnchorShiftBeyondMaxIndelShift_drops() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("GTATCGTGT:NNNN:0:1");

        // Anchor shifted by 2 bases, but maxIndelShift is only 1.
        FastqRecord rec = makeRecord("r1", "XX" + "GTATCGTGT" + "AACC" + "PAYLOAD");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Anchor mismatch tolerance

    @Test
    public void testAnchorOneMismatch_toleratedWithinBudget() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("GTATCGTGT:NNNN:1:0");

        // T -> A at position 4 of the anchor: 1 mismatch.
        FastqRecord rec = makeRecord("r1", "GTATAGTGT" + "AACC" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:AACC", result.getName());
    }

    @Test
    public void testAnchorTwoMismatches_beyondBudget_drops() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("GTATCGTGT:NNNN:1:0");

        // Two substitutions in the anchor -- beyond maxMismatch=1.
        FastqRecord rec = makeRecord("r1", "GAATAGTGT" + "AACC" + "PAYLOAD");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Structured degenerate UMI pattern validation

    @Test
    public void testStructuredPattern_validUmiSurvives() {
        // V = A/C/G only (matches ONT/Karst-style pattern excluding T-homopolymer risk)
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:VVVV:0:0");

        FastqRecord rec = makeRecord("r1", "ACGA" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:ACGA", result.getName());
    }

    @Test
    public void testStructuredPattern_tBaseViolatesVCode_dropsAtZeroMismatch() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:VVVV:0:0");

        // A T at a V-only position is a pattern mismatch.
        FastqRecord rec = makeRecord("r1", "ACGT" + "PAYLOAD");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testStructuredPattern_toleratedWithinMismatchBudget() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:VVVV:1:0");

        FastqRecord rec = makeRecord("r1", "ACGT" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_UMI:ACGT", result.getName());
    }

    // ------------------------------------------------------------------
    // Length edge cases

    @Test
    public void testReadShorterThanUmi_drops() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNNNNNN:0:0");

        FastqRecord rec = makeRecord("r1", "ACGT"); // shorter than 8bp UMI
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testExactLengthRead_zeroLengthPayloadSurvives() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNN:0:0");

        FastqRecord rec = makeRecord("r1", "ACGT");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Custom separator

    @Test
    public void testCustomSeparator() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNN:0:0:__");

        FastqRecord rec = makeRecord("r1", "ACGT" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1__UMI:ACGT", result.getName());
    }

    // ------------------------------------------------------------------
    // Argument validation

    @Test
    public void testMissingRequiredArgs_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("NONE:NNNN"));
    }

    @Test
    public void testEmptyPattern_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("NONE::0:0"));
    }

    @Test
    public void testUnsupportedPatternSymbol_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("NONE:NNNZ:0:0"));
    }

    @Test
    public void testAnchorWithNonBaseChars_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("GTAXT:NNNN:0:0"));
    }

    @Test
    public void testMaxIndelShiftAboveCap_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("NONE:NNNN:0:6"));
    }

    @Test
    public void testNegativeMaxMismatch_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiLongReadExtractTrimmer("NONE:NNNN:-1:0"));
    }

    // ------------------------------------------------------------------
    // Symmetric-mode guard

    @Test
    public void testProcessRecords_twoElementArray_throws() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNN:0:0");

        FastqRecord a = makeRecord("a", "ACGTPAYLOAD");
        FastqRecord b = makeRecord("b", "ACGTPAYLOAD");
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(new FastqRecord[] { a, b }));
    }

    @Test
    public void testProcessRecords_singleElementArray_works() {
        UmiLongReadExtractTrimmer trimmer = new UmiLongReadExtractTrimmer("NONE:NNNN:0:0");

        FastqRecord[] in = { makeRecord("r1", "ACGTPAYLOAD") };
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNotNull(result[0]);
    }
}
