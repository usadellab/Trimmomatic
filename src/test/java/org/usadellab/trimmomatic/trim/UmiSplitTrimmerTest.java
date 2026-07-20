package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class UmiSplitTrimmerTest {

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    @Test
    public void testBasicSplit() {
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("16:12");
        // 16bp CB + 12bp UMI + 8bp payload = 36bp
        String cb = "AAACCCAAGAAACACT";
        String umi = "ACGTACGTACGT";
        String payload = "TTTTGGGG";
        FastqRecord rec = makeRecord("read1", cb + umi + payload);

        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("read1_CB:" + cb + "_UMI:" + umi, result.getName());
        assertEquals(payload, result.getSequence());
    }

    @Test
    public void testCustomSeparator() {
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("4:4:__");
        FastqRecord rec = makeRecord("r1", "AAAACCCCTTTT");

        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1__CB:AAAA__UMI:CCCC", result.getName());
        assertEquals("TTTT", result.getSequence());
    }

    @Test
    public void testExactLengthRead_zeroLengthPayloadSurvives() {
        // Unlike UMIEXTRACT, no leftover payload is required.
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("4:4");
        FastqRecord rec = makeRecord("r1", "AAAACCCC"); // exactly cbLength+umiLength

        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_CB:AAAA_UMI:CCCC", result.getName());
        assertEquals("", result.getSequence());
        assertEquals(0, result.getLength());
    }

    @Test
    public void testShorterThanConstruct_drops() {
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("4:4");
        FastqRecord rec = makeRecord("r1", "AAAACCC"); // 7bp, one short of 8
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testProcessRecords_singleElementArray_works() {
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("4:4");
        FastqRecord[] in = { makeRecord("r1", "AAAACCCCTTTT") };
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNotNull(result[0]);
        assertEquals("TTTT", result[0].getSequence());
    }

    @Test
    public void testProcessRecords_twoElementArray_throws() {
        UmiSplitTrimmer trimmer = new UmiSplitTrimmer("4:4");
        FastqRecord[] in = { makeRecord("a", "AAAACCCCTTTT"), makeRecord("b", "AAAACCCCTTTT") };
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(in));
    }

    @Test
    public void testMissingUmiLength_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiSplitTrimmer("16"));
    }

    @Test
    public void testNegativeCbLength_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiSplitTrimmer("-1:12"));
    }

    @Test
    public void testNegativeUmiLength_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiSplitTrimmer("16:-1"));
    }
}
