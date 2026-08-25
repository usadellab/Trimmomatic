package org.usadellab.trimmomatic.threading.parser;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

/**
 * Tests for InterleavedParserPair.
 *
 * The class reads a single FASTQ file containing alternating R1/R2 records and
 * exposes two virtual Parser views.
 *
 * Covered here:
 *   - Empty input file → immediate EOF on both parsers
 *   - Exactly one pair → both parsers deliver one record each, then EOF
 *   - Multiple pairs across blocks → records ordered and paired correctly
 *   - Odd record count → exception is surfaced (via close() and via poll())
 *   - R1 and R2 sequences are in the correct positions
 *   - Quality strings are preserved per-read
 *   - close() re-throws a stored reader-thread exception even after the thread dies
 */
public class InterleavedParserPairTest {

    @TempDir
    Path tempDir;

    // ------------------------------------------------------------------
    // Helpers

    /** Write an interleaved FASTQ file and return its path. */
    private File writeFastq(String... records) throws Exception {
        File f = tempDir.resolve("interleaved.fastq").toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String record : records) {
                fw.write(record);
            }
        }
        return f;
    }

    /** Build a minimal FASTQ record as a string block (4 lines). */
    private static String fqRecord(String name, String seq) {
        String qual = "I".repeat(seq.length());
        return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
    }

    /** Open a FastqParser with fixed phredOffset=33 to skip the 10k-record pre-read. */
    private FastqParser openParser(File f) throws Exception {
        FastqParser p = new FastqParser(33);
        p.open(f);
        return p;
    }

    /**
     * Poll a parser repeatedly until a non-null list is returned.
     * Honours the 100 ms poll timeout in HalfParser.
     */
    private List<FastqRecord> pollBlock(Parser parser) throws Exception {
        List<FastqRecord> block = null;
        while (block == null) {
            block = parser.poll();
        }
        return block;
    }

    /**
     * Drain all blocks from a parser until the EOF sentinel (empty list).
     * Returns all records received before EOF.
     */
    private List<FastqRecord> drainAll(Parser parser) throws Exception {
        List<FastqRecord> all = new ArrayList<>();
        while (true) {
            List<FastqRecord> block = pollBlock(parser);
            if (block.isEmpty()) {
                break;
            }
            all.addAll(block);
        }
        return all;
    }

    // ------------------------------------------------------------------
    // Empty input file

    @Test
    public void testEmptyFile_bothParsersImmediatelyEof() throws Exception {
        File f = writeFastq(); // empty
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1Block = pollBlock(pair.getR1Parser());
        List<FastqRecord> r2Block = pollBlock(pair.getR2Parser());

        assertTrue(r1Block.isEmpty(), "R1 should be empty (EOF) for empty input");
        assertTrue(r2Block.isEmpty(), "R2 should be empty (EOF) for empty input");

        pair.close();
    }

    // ------------------------------------------------------------------
    // Exactly one pair

    @Test
    public void testOnePair_deliveredCorrectly() throws Exception {
        File f = writeFastq(
                fqRecord("r1/1", "ACGTACGT"),
                fqRecord("r1/2", "TGCATGCA")
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1All = drainAll(pair.getR1Parser());
        List<FastqRecord> r2All = drainAll(pair.getR2Parser());

        assertEquals(1, r1All.size());
        assertEquals(1, r2All.size());
        assertEquals("ACGTACGT", r1All.get(0).getSequence());
        assertEquals("TGCATGCA", r2All.get(0).getSequence());

        pair.close();
    }

    @Test
    public void testOnePair_namesPreserved() throws Exception {
        File f = writeFastq(
                fqRecord("read1/1", "AAAA"),
                fqRecord("read1/2", "CCCC")
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1All = drainAll(pair.getR1Parser());
        List<FastqRecord> r2All = drainAll(pair.getR2Parser());

        assertEquals("read1/1", r1All.get(0).getName());
        assertEquals("read1/2", r2All.get(0).getName());

        pair.close();
    }

    // ------------------------------------------------------------------
    // Multiple pairs, ordering and correct pairing

    @Test
    public void testFourPairs_correctOrderAndPairing() throws Exception {
        File f = writeFastq(
                fqRecord("r1/1", "AAAA"),
                fqRecord("r1/2", "TTTT"),
                fqRecord("r2/1", "CCCC"),
                fqRecord("r2/2", "GGGG"),
                fqRecord("r3/1", "ACGT"),
                fqRecord("r3/2", "TGCA"),
                fqRecord("r4/1", "GGGGA"),
                fqRecord("r4/2", "CCCCT")
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1All = drainAll(pair.getR1Parser());
        List<FastqRecord> r2All = drainAll(pair.getR2Parser());

        assertEquals(4, r1All.size());
        assertEquals(4, r2All.size());

        // R1: odd-indexed records
        assertEquals("AAAA",  r1All.get(0).getSequence());
        assertEquals("CCCC",  r1All.get(1).getSequence());
        assertEquals("ACGT",  r1All.get(2).getSequence());
        assertEquals("GGGGA", r1All.get(3).getSequence());

        // R2: even-indexed records
        assertEquals("TTTT",  r2All.get(0).getSequence());
        assertEquals("GGGG",  r2All.get(1).getSequence());
        assertEquals("TGCA",  r2All.get(2).getSequence());
        assertEquals("CCCCT", r2All.get(3).getSequence());

        pair.close();
    }

    @Test
    public void testBlocksContainEqualNumberOfRecords() throws Exception {
        // 6 pairs; block size is 32768 so they all land in one block.
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 6; i++) {
            sb.append(fqRecord("r" + i + "/1", "ACGT"));
            sb.append(fqRecord("r" + i + "/2", "TGCA"));
        }
        File f = writeFastq(sb.toString());
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        // The first data block from R1 and R2 should have the same size.
        List<FastqRecord> r1Block = pollBlock(pair.getR1Parser());
        List<FastqRecord> r2Block = pollBlock(pair.getR2Parser());

        assertEquals(r1Block.size(), r2Block.size(),
                "R1 and R2 data blocks must have equal sizes");

        pair.close();
    }

    @Test
    public void testQualityStringsPreservedPerRead() throws Exception {
        // Use distinct quality strings so we can tell R1 from R2.
        File f = tempDir.resolve("qual_test.fastq").toFile();
        try (FileWriter fw = new FileWriter(f)) {
            fw.write("@r1/1\nACGT\n+\nIIII\n");
            fw.write("@r1/2\nTGCA\n+\nJJJJ\n");
        }
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1All = drainAll(pair.getR1Parser());
        List<FastqRecord> r2All = drainAll(pair.getR2Parser());

        assertEquals("IIII", r1All.get(0).getQuality(), "R1 quality should be IIII");
        assertEquals("JJJJ", r2All.get(0).getQuality(), "R2 quality should be JJJJ");

        pair.close();
    }

    // ------------------------------------------------------------------
    // Odd record count → exception

    @Test
    public void testOddRecordCount_exceptionSurfacedOnClose() throws Exception {
        // 3 records: pair1 + orphan R1 record, must throw
        File f = writeFastq(
                fqRecord("r1/1", "AAAA"),
                fqRecord("r1/2", "CCCC"),
                fqRecord("r2/1", "GGGG")  // no partner
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        // Drain everything the parsers offer (may include data from pair1)
        drainAllIgnoringExceptions(pair.getR1Parser());
        drainAllIgnoringExceptions(pair.getR2Parser());

        // close() must rethrow, even if the reader thread has already died
        assertThrows(Exception.class, pair::close);
    }

    @Test
    public void testOddRecordCount_exceptionSurfacedOnPoll() throws Exception {
        File f = writeFastq(
                fqRecord("r1/1", "AAAA"),
                fqRecord("r1/2", "CCCC"),
                fqRecord("r2/1", "GGGG")
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        // Poll until we either hit the exception or exhaust the data.
        boolean exceptionSeen = false;
        outer:
        for (int round = 0; round < 200; round++) {
            for (Parser parser : new Parser[]{ pair.getR1Parser(), pair.getR2Parser() }) {
                try {
                    List<FastqRecord> block = parser.poll();
                    if (block != null && block.isEmpty()) {
                        // Received EOF sentinel without exception, reader may not have
                        // failed yet; keep retrying for a short period.
                        Thread.sleep(5);
                    }
                } catch (Exception e) {
                    exceptionSeen = true;
                    break outer;
                }
            }
        }

        assertTrue(exceptionSeen, "Expected exception to surface during poll for odd-record file");

        // close() will re-throw, but source.close() is guaranteed to run (try-finally).
        try { pair.close(); } catch (Exception ignored) {}
    }

    @Test
    public void testSingleRecord_exceptionSurfacedOnClose() throws Exception {
        // 1 record (minimum possible odd count)
        File f = writeFastq(fqRecord("r1/1", "ACGT"));
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        drainAllIgnoringExceptions(pair.getR1Parser());
        drainAllIgnoringExceptions(pair.getR2Parser());

        assertThrows(Exception.class, pair::close);
    }

    // ------------------------------------------------------------------
    // Different read lengths in the same pair

    @Test
    public void testMixedLengthReads_deliveredCorrectly() throws Exception {
        File f = writeFastq(
                fqRecord("r1/1", "ACGT"),                           // 4 bases
                fqRecord("r1/2", "TGCATGCATGCATGCATGCATGCATGCA")   // 28 bases
        );
        FastqParser src = openParser(f);
        ExceptionHolder eh = new ExceptionHolder();
        InterleavedParserPair pair = new InterleavedParserPair(src, 4, eh);

        List<FastqRecord> r1All = drainAll(pair.getR1Parser());
        List<FastqRecord> r2All = drainAll(pair.getR2Parser());

        assertEquals(4,  r1All.get(0).getLength());
        assertEquals(28, r2All.get(0).getLength());

        pair.close();
    }

    // ------------------------------------------------------------------
    // Private helper: drain ignoring exceptions (used before asserting on close)

    private void drainAllIgnoringExceptions(Parser parser) {
        // Give the reader thread time to finish, then drain whatever was queued.
        try { Thread.sleep(50); } catch (InterruptedException ignored) {}
        for (int i = 0; i < 50; i++) {
            try {
                List<FastqRecord> block = parser.poll();
                if (block != null && block.isEmpty()) return; // EOF sentinel
            } catch (Exception ignored) {
                return; // exception observed, stop draining
            }
        }
    }
}
