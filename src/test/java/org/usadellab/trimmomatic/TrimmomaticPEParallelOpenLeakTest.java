package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

/**
 * TrimmomaticPE.process() opens both R1 and R2 in parallel on virtual threads.
 * If R2's open() fails after R1's already succeeded, R1's FastqParser used to
 * never get closed. On Windows a file with a live read handle cannot be
 * deleted, so that leak is directly observable.
 */
public class TrimmomaticPEParallelOpenLeakTest {

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
    public void testSuccessfullyOpenedParserIsClosedWhenSiblingOpenFails(@TempDir File tempDir) throws IOException {
        File goodInput = writeFastq(tempDir, "good.fastq", fqRecord("r1", "ACGT"));
        File badInput = new File(tempDir, "does-not-exist.fastq");

        File r1p = new File(tempDir, "r1p.fastq");
        File r1u = new File(tempDir, "r1u.fastq");
        File r2p = new File(tempDir, "r2p.fastq");
        File r2u = new File(tempDir, "r2u.fastq");

        Logger logger = new Logger(false, false, false);
        TrimmomaticPE tm = new TrimmomaticPE(logger);

        assertThrows(Exception.class, () -> tm.process(goodInput, badInput, false, r1p, r1u, r2p, r2u,
                new Trimmer[0], 33, null, null, false, null, null, 1, false, 0));

        assertTrue(goodInput.delete(),
                "Successfully-opened input file must be closed after the sibling's open() failure");
    }
}
