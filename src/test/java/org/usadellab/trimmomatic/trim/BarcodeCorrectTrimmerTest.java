package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BarcodeCorrectTrimmerTest {

    @TempDir
    Path tempDir;

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    private File writeWhitelist(String filename, String... barcodes) throws IOException {
        File f = tempDir.resolve(filename).toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String b : barcodes) fw.write(b + "\n");
        }
        return f;
    }

    private File writeGzippedWhitelist(String filename, String... barcodes) throws IOException {
        File f = tempDir.resolve(filename).toFile();
        try (GZIPOutputStream gz = new GZIPOutputStream(new java.io.FileOutputStream(f))) {
            for (String b : barcodes)
                gz.write((b + "\n").getBytes());
        }
        return f;
    }

    // ------------------------------------------------------------------
    // Exact match

    @Test
    public void testExactMatch_survivesWithSameBarcode() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord rec = makeRecord("r1", "AAAACCCCGGGGTTTT" + "UMIUMIUM" + "PAYLOAD1");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AAAACCCCGGGGTTTT_UMIUMIUM", result.getName());
        assertEquals("PAYLOAD1", result.getSequence());
    }

    // ------------------------------------------------------------------
    // 1-mismatch correction

    @Test
    public void testOneMismatch_correctsToWhitelistEntry() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        // Position 0: A -> T (1 mismatch from the whitelist entry)
        String rawCb = "TAAACCCCGGGGTTTT";
        FastqRecord rec = makeRecord("r1", rawCb + "UMIUMIUM" + "PAYLOAD1");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AAAACCCCGGGGTTTT_UMIUMIUM", result.getName(),
                "corrected CB in the name must be the WHITELIST entry, not the raw sequence");
    }

    @Test
    public void testNInBarcode_correctable() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        String rawCb = "NAAACCCCGGGGTTTT"; // N at position 0
        FastqRecord rec = makeRecord("r1", rawCb + "UMIUMIUM" + "PAYLOAD1");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AAAACCCCGGGGTTTT_UMIUMIUM", result.getName());
    }

    // ------------------------------------------------------------------
    // No match within tolerance

    @Test
    public void testTwoMismatches_beyondTolerance_drops() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        // Positions 0 and 1 both wrong -- 2 mismatches, beyond maxMismatch=1
        String rawCb = "TTAACCCCGGGGTTTT";
        FastqRecord rec = makeRecord("r1", rawCb + "UMIUMIUM" + "PAYLOAD1");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testMaxMismatchZero_onlyExactMatchWorks() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:0");

        String oneOff = "TAAACCCCGGGGTTTT"; // 1 mismatch -- should NOT correct at maxMismatch=0
        FastqRecord rec = makeRecord("r1", oneOff + "UMIUMIUM" + "PAYLOAD1");
        assertNull(trimmer.processRecord(rec));

        String exact = "AAAACCCCGGGGTTTT";
        FastqRecord rec2 = makeRecord("r2", exact + "UMIUMIUM" + "PAYLOAD1");
        assertNotNull(trimmer.processRecord(rec2));
    }

    // ------------------------------------------------------------------
    // Length edge cases

    @Test
    public void testExactLengthRead_zeroLengthPayloadSurvives() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord rec = makeRecord("r1", "AAAACCCCGGGGTTTT" + "UMIUMIUM"); // no payload
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("", result.getSequence());
        assertEquals(0, result.getLength());
    }

    @Test
    public void testShorterThanConstruct_drops() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord rec = makeRecord("r1", "AAAACCCCGGGGTTT"); // 15bp, one short of 16+8
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Whitelist file loading

    @Test
    public void testGzippedWhitelist_loadsCorrectly() throws Exception {
        File wl = writeGzippedWhitelist("wl.txt.gz", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord rec = makeRecord("r1", "AAAACCCCGGGGTTTT" + "UMIUMIUM" + "PAYLOAD1");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testEmptyWhitelist_throws() throws Exception {
        File wl = writeWhitelist("empty.txt");
        assertThrows(IllegalArgumentException.class,
                () -> new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1"));
    }

    // ------------------------------------------------------------------
    // Argument validation

    @Test
    public void testMismatchAboveSupportedCap_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        assertThrows(IllegalArgumentException.class,
                () -> new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:2"));
    }

    @Test
    public void testMissingRequiredArgs_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        assertThrows(IllegalArgumentException.class,
                () -> new BarcodeCorrectTrimmer(wl.getPath() + ":16"));
    }

    // ------------------------------------------------------------------
    // Symmetric-mode guard

    @Test
    public void testProcessRecords_twoElementArray_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord a = makeRecord("a", "AAAACCCCGGGGTTTT" + "UMIUMIUM" + "X");
        FastqRecord b = makeRecord("b", "AAAACCCCGGGGTTTT" + "UMIUMIUM" + "X");
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(new FastqRecord[] { a, b }));
    }

    @Test
    public void testProcessRecords_singleElementArray_works() throws Exception {
        File wl = writeWhitelist("wl.txt", "AAAACCCCGGGGTTTT");
        BarcodeCorrectTrimmer trimmer = new BarcodeCorrectTrimmer(wl.getPath() + ":16:8:1");

        FastqRecord[] in = { makeRecord("r1", "AAAACCCCGGGGTTTT" + "UMIUMIUM" + "PAYLOAD1") };
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNotNull(result[0]);
    }
}
