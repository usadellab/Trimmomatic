package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class UmiDimerCorrectTrimmerTest {

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

    // ------------------------------------------------------------------
    // Exact match

    @Test
    public void testExactMatch_survives() throws Exception {
        // 4 dimers (8bp) cell barcode: AA-CC-GG-TT
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        FastqRecord rec = makeRecord("r1", "AACCGGTT" + "UMIU" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AACCGGTT_UMIU", result.getName());
        assertEquals("PAYLOAD", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Single-dimer correction absorbing a 2-base error

    @Test
    public void testSingleDimerBothBasesWrong_corrected() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        // First dimer AA -> TT (both bases wrong), rest untouched, one corrupted dimer.
        String rawCb = "TTCCGGTT";
        FastqRecord rec = makeRecord("r1", rawCb + "UMIU" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AACCGGTT_UMIU", result.getName(),
                "corrected CB must be the WHITELIST entry, not the raw sequence");
    }

    @Test
    public void testSingleDimerOneBaseWrong_corrected() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        // Third dimer GG -> GA (one base wrong).
        String rawCb = "AACCGATT";
        FastqRecord rec = makeRecord("r1", rawCb + "UMIU" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AACCGGTT_UMIU", result.getName());
    }

    // ------------------------------------------------------------------
    // Two corrupted dimers, beyond correction capability

    @Test
    public void testTwoDimersCorrupted_drops() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        // Dimer 1 (AA->TT) AND dimer 3 (GG->CC) both corrupted.
        String rawCb = "TTCCCCTT";
        FastqRecord rec = makeRecord("r1", rawCb + "UMIU" + "PAYLOAD");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testMaxMismatchDimersZero_onlyExactMatchWorks() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:0");

        FastqRecord oneOff = makeRecord("r1", "TTCCGGTT" + "UMIU" + "PAYLOAD");
        assertNull(trimmer.processRecord(oneOff));

        FastqRecord exact = makeRecord("r2", "AACCGGTT" + "UMIU" + "PAYLOAD");
        assertNotNull(trimmer.processRecord(exact));
    }

    // ------------------------------------------------------------------
    // UMI never corrected, passed through raw

    @Test
    public void testUmiPassedThroughRaw_notWhitelistChecked() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        FastqRecord rec = makeRecord("r1", "AACCGGTT" + "ZZZZ" + "PAYLOAD");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("r1_AACCGGTT_ZZZZ", result.getName());
    }

    // ------------------------------------------------------------------
    // Length edge cases

    @Test
    public void testExactLengthRead_zeroLengthPayloadSurvives() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        FastqRecord rec = makeRecord("r1", "AACCGGTT" + "UMIU");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("", result.getSequence());
    }

    @Test
    public void testShorterThanConstruct_drops() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        FastqRecord rec = makeRecord("r1", "AACCGGTTUMI"); // 11bp, one short of 8+4
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Argument validation

    @Test
    public void testOddCbLength_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        assertThrows(IllegalArgumentException.class,
                () -> new UmiDimerCorrectTrimmer(wl.getPath() + ":7:4:1"));
    }

    @Test
    public void testMismatchDimersAboveSupportedCap_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        assertThrows(IllegalArgumentException.class,
                () -> new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:2"));
    }

    @Test
    public void testEmptyWhitelist_throws() throws Exception {
        File wl = writeWhitelist("empty.txt");
        assertThrows(IllegalArgumentException.class,
                () -> new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1"));
    }

    @Test
    public void testMissingRequiredArgs_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        assertThrows(IllegalArgumentException.class,
                () -> new UmiDimerCorrectTrimmer(wl.getPath() + ":8"));
    }

    // ------------------------------------------------------------------
    // Symmetric-mode guard

    @Test
    public void testProcessRecords_twoElementArray_throws() throws Exception {
        File wl = writeWhitelist("wl.txt", "AACCGGTT");
        UmiDimerCorrectTrimmer trimmer = new UmiDimerCorrectTrimmer(wl.getPath() + ":8:4:1");

        FastqRecord a = makeRecord("a", "AACCGGTT" + "UMIU" + "X");
        FastqRecord b = makeRecord("b", "AACCGGTT" + "UMIU" + "X");
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(new FastqRecord[] { a, b }));
    }
}
