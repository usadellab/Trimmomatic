package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Tests for UmiExtractTrimmer.processRecords(FastqRecord[]), the entry
 * point BlockOfWork actually calls, as opposed to processRecord(FastqRecord)
 * which the other UmiExtractTrimmer test classes exercise directly.
 *
 * UMIEXTRACT renames a record from its OWN leading bases. Handed a 2-element
 * array (symmetric PE mode: both mates in one call), it would previously
 * rename both mates independently from their own, different leading bases.
 * That desynchronised mate names with no error raised, and chopped real bases
 * off whichever mate is not the barcode/UMI read. It now refuses that call
 * shape outright instead of corrupting data.
 */
public class UmiExtractTrimmerProcessRecordsTest {

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    @Test
    public void testSingleElementArray_worksNormally() {
        // This is the correct call shape: one mate, in its own single-element
        // array, as used by -pe1steps/-pe2steps and -technicalread.
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord[] in = { makeRecord("read1", "TTTGACGT") };
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNotNull(result[0]);
        assertEquals("read1_UMI:TTTG", result[0].getName());
        assertEquals("ACGT", result[0].getSequence());
    }

    @Test
    public void testSingleElementArray_dropReturnsNullElement() {
        // Regression: the single-element path must still preserve the existing
        // drop-on-too-short behaviour (array element becomes null, not an
        // exception, not a shorter array).
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("8");
        FastqRecord[] in = { makeRecord("read1", "ACGTACGT") }; // length == umiLength -> drop
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNull(result[0]);
    }

    @Test
    public void testNullArray_returnsNull() {
        // AbstractSingleRecordTrimmer.processRecords(null) returns null; the
        // guard must not NPE on the length check before delegating.
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        assertNull(trimmer.processRecords(null));
    }

    @Test
    public void testTwoElementArray_symmetricPEMode_throws() {
        // This is exactly the call shape BlockOfWork's symmetric (non-per-mate,
        // non-technicalRead) PE branch uses: both mates in one array. Must
        // refuse, instead of renaming and corrupting both independently.
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord[] in = {
                makeRecord("pair/1", "TTTGACGTACGT"),
                makeRecord("pair/2", "GGGGCCCCTTTT")
        };

        IllegalStateException ex = assertThrows(IllegalStateException.class,
                () -> trimmer.processRecords(in));
        // Message should point users at the actual fix, not just say "no".
        assertEquals(true, ex.getMessage().contains("-pe1steps"));
        assertEquals(true, ex.getMessage().contains("-technicalread"));
    }

    @Test
    public void testThreeElementArray_alsoThrows() {
        // Not a real Trimmomatic call shape, but the guard is a plain length
        // check ("> 1"), not a hardcoded "== 2", confirm it generalises.
        UmiExtractTrimmer trimmer = new UmiExtractTrimmer("4");
        FastqRecord[] in = {
                makeRecord("a", "TTTGACGT"),
                makeRecord("b", "TTTGACGT"),
                makeRecord("c", "TTTGACGT")
        };
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(in));
    }
}
