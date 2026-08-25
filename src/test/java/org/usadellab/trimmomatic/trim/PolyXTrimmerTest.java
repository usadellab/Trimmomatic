package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class PolyXTrimmerTest {

    private FastqRecord makeRecord(String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    @Test
    public void testShortPolyATailNotTrimmed() {
        // 5-base polyA tail, minLength=10 → pass unchanged
        PolyXTrimmer trimmer = new PolyXTrimmer("A:10");
        FastqRecord rec = makeRecord("ACGTACGTAAAAAA");
        // tail has 6 A's, still below min 10 → unchanged
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals(rec.getSequence(), result.getSequence());
    }

    @Test
    public void testLongPolyATailTrimmed() {
        // "ACGTACGT" + 12 A's, minLength=10
        PolyXTrimmer trimmer = new PolyXTrimmer("A:10");
        FastqRecord rec = makeRecord("ACGTACGTAAAAAAAAAAAA");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testExactMinLengthTrimmed() {
        // Exactly minLength A's
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("GCGCAAAAA");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    @Test
    public void testEntireReadPolyXDropped() {
        // Whole read is poly-A → drop
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("AAAAAAAAAA");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testEmptyReadDropped() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:5");
        FastqRecord rec = makeRecord("");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testPolyGTrimming() {
        // Common two-colour Illumina artefact: polyG at 3' end
        PolyXTrimmer trimmer = new PolyXTrimmer("G:10");
        FastqRecord rec = makeRecord("ACGTACGTGGGGGGGGGGGG");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testPolyTTrimming() {
        PolyXTrimmer trimmer = new PolyXTrimmer("T:8");
        FastqRecord rec = makeRecord("AAACCCTTTTTTTTT");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("AAACCC", result.getSequence());
    }

    @Test
    public void testNoTailAtAllPassesThrough() {
        PolyXTrimmer trimmer = new PolyXTrimmer("A:10");
        FastqRecord rec = makeRecord("GCGCGCGCGCGCGCGC");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals(rec.getSequence(), result.getSequence());
    }

    @Test
    public void testMissingSecondArgThrows() {
        assertThrows(IllegalArgumentException.class, () -> new PolyXTrimmer("A"));
    }

    @Test
    public void testInvalidBaseThrows() {
        assertThrows(IllegalArgumentException.class, () -> new PolyXTrimmer("X:10"));
    }

    @Test
    public void testZeroMinLengthThrows() {
        assertThrows(IllegalArgumentException.class, () -> new PolyXTrimmer("A:0"));
    }

    @Test
    public void testBaseNAllowed() {
        // N is a valid target base per the spec
        PolyXTrimmer trimmer = new PolyXTrimmer("N:5");
        FastqRecord rec = makeRecord("ACGTNNNNNNNN");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGT", result.getSequence());
    }

    @Test
    public void testLowercaseBaseNormalized() {
        // Parser should upper-case the base spec
        PolyXTrimmer trimmer = new PolyXTrimmer("a:5");
        FastqRecord rec = makeRecord("GCGCAAAAAA");
        FastqRecord result = trimmer.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }
}
