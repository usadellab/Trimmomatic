package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

/**
 * Regression test for a stats-aggregation race in the multi-threaded PE path.
 *
 * SelfThreadedTrimStatsCollector merges each block's TrimStats on its own
 * background thread. TrimmomaticPE.processPipeline() used to read the
 * aggregated stats object (statsCollector.getStats()) immediately after the
 * main reading loop finished submitting all blocks, without first waiting
 * for that background thread to actually finish draining and merging them.
 * For a fast-to-read input the main thread could win that race and log an
 * undercount, confirmed empirically: a 10,000-pair input reproducibly
 * reported "Input Read Pairs: 0" at -threads 4 despite the actual FASTQ
 * output being complete and correct (this bug never affected output
 * correctness, only the printed/written summary numbers).
 *
 * Only reproduces through SelfThreadedTrimStatsCollector, which is only used
 * when threads > 1, ParasiteTrimStatsCollector (threads == 1) merges
 * synchronously on the caller's thread and never raced.
 *
 * Fix: TrimmomaticPE now calls statsCollector.close() (which waits for the
 * background thread to finish) before reading .getStats().
 */
public class TrimmomaticPEThreadedStatsTest {

    @TempDir
    Path tempDir;

    private static final int READ_COUNT = 10000;

    private File writeFastq(String filename, int count) throws Exception {
        File f = tempDir.resolve(filename).toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (int i = 0; i < count; i++) {
                fw.write("@read" + i + "\n");
                fw.write("ACGTACGTACGTACGTACGTACGTACGT\n");
                fw.write("+\n");
                fw.write("IIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
            }
        }
        return f;
    }

    @Test
    public void testStatsCorrectAtThreads4() throws Exception {
        File r1 = writeFastq("r1.fastq", READ_COUNT);
        File r2 = writeFastq("r2.fastq", READ_COUNT);
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();
        File summary = tempDir.resolve("summary.txt").toFile();

        Logger logger = new Logger(false, false, false);
        TrimmomaticPE tm = new TrimmomaticPE(logger);
        tm.process(r1, r2, false, r1p, r1u, r2p, r2u, new Trimmer[0],
                33, null, summary, false, null, null, 4, false, 0);

        List<String> lines = Files.readAllLines(summary.toPath());
        String inputLine = lines.stream().filter(l -> l.startsWith("Input Read Pairs:")).findFirst().orElse(null);

        assertNotNull(inputLine, "summary file must contain an Input Read Pairs line");
        assertEquals("Input Read Pairs: " + READ_COUNT, inputLine,
                "stats must reflect all " + READ_COUNT + " input pairs, not a partial async snapshot");
    }

    @Test
    public void testStatsCorrectAtThreads8() throws Exception {
        // Same check at a different thread count, the race is timing-dependent,
        // not tied to one specific thread count.
        File r1 = writeFastq("r1_t8.fastq", READ_COUNT);
        File r2 = writeFastq("r2_t8.fastq", READ_COUNT);
        File r1p = tempDir.resolve("r1p_t8.fastq").toFile();
        File r1u = tempDir.resolve("r1u_t8.fastq").toFile();
        File r2p = tempDir.resolve("r2p_t8.fastq").toFile();
        File r2u = tempDir.resolve("r2u_t8.fastq").toFile();
        File summary = tempDir.resolve("summary_t8.txt").toFile();

        Logger logger = new Logger(false, false, false);
        TrimmomaticPE tm = new TrimmomaticPE(logger);
        tm.process(r1, r2, false, r1p, r1u, r2p, r2u, new Trimmer[0],
                33, null, summary, false, null, null, 8, false, 0);

        List<String> lines = Files.readAllLines(summary.toPath());
        String inputLine = lines.stream().filter(l -> l.startsWith("Input Read Pairs:")).findFirst().orElse(null);

        assertNotNull(inputLine, "summary file must contain an Input Read Pairs line");
        assertEquals("Input Read Pairs: " + READ_COUNT, inputLine,
                "stats must reflect all " + READ_COUNT + " input pairs, not a partial async snapshot");
    }
}
