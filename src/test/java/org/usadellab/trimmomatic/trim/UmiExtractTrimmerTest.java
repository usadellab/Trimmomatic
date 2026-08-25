package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class UmiExtractTrimmerTest {

    private FastqRecord makeRecord(String name, String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord(name, seq, "", qual, 33);
    }

    // ------------------------------------------------------------------
    // Name modification

    @Test
    public void testUmiAppendedToName() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("6");
        FastqRecord rec = makeRecord("read1", "ACGTACGTACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("read1_UMI:ACGTAC", result.getName());
    }

    @Test
    public void testCustomSeparator() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4:__");
        FastqRecord rec = makeRecord("read1", "TTTGACGT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("read1__UMI:TTTG", result.getName());
    }

    // ------------------------------------------------------------------
    // Sequence after extraction

    @Test
    public void testSequenceTrimmedCorrectly() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("6");
        FastqRecord rec = makeRecord("read1", "UUUUUUACGTACGT");
        // First 6 bases extracted as UMI, remaining = "ACGTACGT"
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testQualityAlsoTrimmedCorrectly() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        // Construct record with distinguishable quality
        String seq = "ACGTACGTACGT";
        String qual = "IIII" + "JJJJJJJJ"; // first 4 = UMI quality, rest = payload quality
        FastqRecord rec = new FastqRecord("read1", seq, "", qual, 33);
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("JJJJJJJJ", result.getQuality());
    }

    // ------------------------------------------------------------------
    // Length checks

    @Test
    public void testReadTooShortDropped() {
        // Read length equals UMI length, no payload remains → drop
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("8");
        FastqRecord rec = makeRecord("read1", "ACGTACGT");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testReadShorterThanUmiDropped() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("12");
        FastqRecord rec = makeRecord("read1", "ACGT");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testReadOneBaseMoreThanUmiKept() {
        // One base of payload remains
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("5");
        FastqRecord rec = makeRecord("read1", "NNNNNX");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("X", result.getSequence());
        assertEquals(1, result.getLength());
    }

    // ------------------------------------------------------------------
    // Lengths and position preservation

    @Test
    public void testResultLengthCorrect() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("10");
        FastqRecord rec = makeRecord("read1", "ACGTACGTACGTACGT"); // 16 bases
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals(6, result.getLength()); // 16 - 10 = 6
    }

    @Test
    public void testUmiIsTakenFrom5Prime() {
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord rec = makeRecord("read1", "TTTTGGGGCCCC");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        // UMI = "TTTT" (first 4 bases), name contains it
        assertTrue(result.getName().contains("TTTT"),
                "Expected UMI 'TTTT' in name, got: " + result.getName());
        // Remaining sequence starts at position 4
        assertEquals("GGGGCCCC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Invalid args

    @Test
    public void testZeroLengthThrows() {
        assertThrows(IllegalArgumentException.class, () -> new UmiExtractTrimmer("0"));
    }

    @Test
    public void testNegativeLengthThrows() {
        assertThrows(IllegalArgumentException.class, () -> new UmiExtractTrimmer("-1"));
    }
}
