package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * Pairomatic.getFastqNames()/splitFastq() used to open a FastqParser and never
 * close it. On Windows a file with a live read handle cannot be deleted, so
 * that leak is directly observable: delete the input file right after
 * process() returns and it must succeed.
 */
public class PairomaticResourceLeakTest {

    private static String fqRecord(String name, String seq) {
        String qual = "I".repeat(seq.length());
        return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
    }

    private File writeFastq(File dir, String filename, String... records) throws IOException {
        File f = new File(dir, filename);
        try (FileWriter fw = new FileWriter(f)) {
            for (String r : records)
                fw.write(r);
        }
        return f;
    }

    @Test
    public void testInputFilesAreClosedAfterProcess(@TempDir File tempDir) throws IOException {
        File in1 = writeFastq(tempDir, "in1.fastq", fqRecord("r1", "ACGT"));
        File in2 = writeFastq(tempDir, "in2.fastq", fqRecord("r1", "ACGT"));
        File out1P = new File(tempDir, "out1P.fastq");
        File out1U = new File(tempDir, "out1U.fastq");
        File out2P = new File(tempDir, "out2P.fastq");
        File out2U = new File(tempDir, "out2U.fastq");

        new Pairomatic().process(in1, in2, out1P, out1U, out2P, out2U, null);

        assertTrue(in1.delete(), "Input file 1 must be closed (and therefore deletable) after process()");
        assertTrue(in2.delete(), "Input file 2 must be closed (and therefore deletable) after process()");
    }
}
