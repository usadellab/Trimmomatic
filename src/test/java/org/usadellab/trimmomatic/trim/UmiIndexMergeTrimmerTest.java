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

public class UmiIndexMergeTrimmerTest {

    @TempDir
    Path tempDir;

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    private File writeIndexFastq(String filename, String... nameSeqPairs) throws IOException {
        File f = tempDir.resolve(filename).toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (int i = 0; i < nameSeqPairs.length; i += 2) {
                String name = nameSeqPairs[i];
                String seq = nameSeqPairs[i + 1];
                fw.write("@" + name + "\n" + seq + "\n+\n" + "I".repeat(seq.length()) + "\n");
            }
        }
        return f;
    }

    // ------------------------------------------------------------------
    // Basic merge by shared read ID

    @Test
    public void testMatchingReadId_tagsUmiFromIndexRead() throws Exception {
        File idx = writeIndexFastq("I2.fastq", "read1", "ACGTACGT");
        UmiIndexMergeTrimmer trimmer = new UmiIndexMergeTrimmer(idx.getPath());

        FastqRecord rec = makeRecord("read1", "PAYLOADSEQ");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("read1_UMI:ACGTACGT", result.getName());
        assertEquals("PAYLOADSEQ", result.getSequence(), "read sequence itself must be untouched, UMI lives in a separate file");
    }

    // ------------------------------------------------------------------
    // Read ID matching ignores the space-delimited description

    @Test
    public void testReadIdMatchIgnoresDescriptionAfterSpace() throws Exception {
        File idx = writeIndexFastq("I2.fastq", "read1 2:N:0:INDEX", "ACGTACGT");
        UmiIndexMergeTrimmer trimmer = new UmiIndexMergeTrimmer(idx.getPath());

        FastqRecord rec = makeRecord("read1 1:N:0:INDEX", "PAYLOADSEQ");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("read1 1:N:0:INDEX_UMI:ACGTACGT", result.getName());
    }

    // ------------------------------------------------------------------
    // Missing index entry

    @Test
    public void testNoMatchingIndexEntry_drops() throws Exception {
        File idx = writeIndexFastq("I2.fastq", "read1", "ACGTACGT");
        UmiIndexMergeTrimmer trimmer = new UmiIndexMergeTrimmer(idx.getPath());

        FastqRecord rec = makeRecord("read2", "PAYLOADSEQ");
        assertNull(trimmer.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Symmetric PE safety: both mates share read ID, get the same UMI tag

    @Test
    public void testSymmetricPairedEnd_bothMatesTaggedIdentically() throws Exception {
        File idx = writeIndexFastq("I2.fastq", "read1", "ACGTACGT");
        UmiIndexMergeTrimmer trimmer = new UmiIndexMergeTrimmer(idx.getPath());

        FastqRecord r1 = makeRecord("read1", "R1SEQ");
        FastqRecord r2 = makeRecord("read1", "R2SEQ");

        FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { r1, r2 });

        assertEquals(2, result.length);
        assertEquals("read1_UMI:ACGTACGT", result[0].getName());
        assertEquals("read1_UMI:ACGTACGT", result[1].getName());
    }

    // ------------------------------------------------------------------
    // Custom separator

    @Test
    public void testCustomSeparator() throws Exception {
        File idx = writeIndexFastq("I2.fastq", "read1", "ACGTACGT");
        UmiIndexMergeTrimmer trimmer = new UmiIndexMergeTrimmer(idx.getPath() + ":__");

        FastqRecord rec = makeRecord("read1", "PAYLOADSEQ");
        FastqRecord result = trimmer.processRecord(rec);

        assertNotNull(result);
        assertEquals("read1__UMI:ACGTACGT", result.getName());
    }

    // ------------------------------------------------------------------
    // Argument validation

    @Test
    public void testEmptyIndexFastq_throws() throws Exception {
        File idx = tempDir.resolve("empty.fastq").toFile();
        idx.createNewFile();
        assertThrows(IllegalArgumentException.class, () -> new UmiIndexMergeTrimmer(idx.getPath()));
    }

    @Test
    public void testMissingArgs_throws() {
        assertThrows(IllegalArgumentException.class, () -> new UmiIndexMergeTrimmer(""));
    }
}
