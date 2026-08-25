package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.trim.HeadCropTrimmer;
import org.usadellab.trimmomatic.trim.MinLenTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

/**
 * Tests for the -technicalread flag in TrimmomaticPE.
 *
 * -technicalread 1 → R1 is the technical read (barcode/UMI), R2 is trimmed.
 * -technicalread 2 → R2 is the technical read, R1 is trimmed.
 *
 * Key semantics verified here:
 *   - The technical read passes through completely untouched.
 *   - All trimmers run only on the biological read.
 *   - If the biological read is dropped, BOTH reads are discarded.
 *   - The technical read NEVER appears in the unpaired output.
 *   - Normal PE behaviour (technicalRead=0) is unaffected.
 */
public class TrimmomaticPETechnicalReadTest {

    @TempDir
    Path tempDir;

    // -----------------------------------------------------------------------
    // Helpers

    private static String fqRecord(String name, String seq) {
        String qual = "I".repeat(seq.length());
        return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
    }

    private File writeFastq(String filename, String... records) throws Exception {
        File f = tempDir.resolve(filename).toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String r : records) fw.write(r);
        }
        return f;
    }

    /** Read all sequence lines (every 4th line starting from line 2) from a FASTQ file. */
    private List<String> readSequences(File f) throws Exception {
        List<String> seqs = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            int lineNo = 0;
            while ((line = br.readLine()) != null) {
                if (lineNo % 4 == 1) seqs.add(line);
                lineNo++;
            }
        }
        return seqs;
    }

    /** Run TrimmomaticPE.process() with the given technicalRead value and trimmers. */
    private void run(File r1in, File r2in, File r1p, File r1u, File r2p, File r2u,
                     int technicalRead, Trimmer... trimmers) throws Exception {
        Logger logger = new Logger(false, false, false);
        TrimmomaticPE tm = new TrimmomaticPE(logger);
        tm.process(r1in, r2in, false, r1p, r1u, r2p, r2u, trimmers,
                33, null, null, false, null, null, 1, false, technicalRead);
    }

    // -----------------------------------------------------------------------
    // technicalRead=1 (R1=tech, R2=bio), bio read survives

    @Test
    public void testTechRead1_bioSurvives_bothInPairedOutput() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGT"));          // 8 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGTACGTACGTACGT")); // 16 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(5));

        assertEquals(List.of("ACGTACGT"),          readSequences(r1p));
        assertEquals(List.of("ACGTACGTACGTACGT"),  readSequences(r2p));
        assertTrue(readSequences(r1u).isEmpty(), "R1 unpaired must be empty");
        assertTrue(readSequences(r2u).isEmpty(), "R2 unpaired must be empty");
    }

    @Test
    public void testTechRead1_techReadSequenceIsUnchanged() throws Exception {
        // Tech read has N's at the start, HeadCrop(4) would strip them if applied.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "NNNNACGT"));    // 8 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "NNNNTTTTTTTT")); // 12 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new HeadCropTrimmer("4"));

        // R1 (tech): must be completely unchanged, 4 N's must still be there
        assertEquals(List.of("NNNNACGT"), readSequences(r1p),
                "Technical read must not be modified by HeadCrop");
        // R2 (bio): first 4 bases cropped
        assertEquals(List.of("TTTTTTTT"), readSequences(r2p));
    }

    @Test
    public void testTechRead1_bioTrimmedBySurvives_techUnchanged() throws Exception {
        // HeadCrop(5): bio goes from 13bp -> 8bp (survives); tech must remain 9bp.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "NNNNNACGT"));        // 9 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "NNNNNTTTTTTTT"));    // 13 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new HeadCropTrimmer("5"));

        assertEquals(List.of("NNNNNACGT"), readSequences(r1p), "Tech read must be 9bp unchanged");
        assertEquals(List.of("TTTTTTTT"),  readSequences(r2p), "Bio read must be trimmed to 8bp");
        assertTrue(readSequences(r1u).isEmpty());
        assertTrue(readSequences(r2u).isEmpty());
    }

    // -----------------------------------------------------------------------
    // technicalRead=1, bio read dropped

    @Test
    public void testTechRead1_bioDropped_allOutputsEmpty() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGTACGTACGT")); // 16 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGT"));              // 4 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(10)); // 4 < 10 → bio dropped

        assertTrue(readSequences(r1p).isEmpty(), "R1 paired must be empty");
        assertTrue(readSequences(r2p).isEmpty(), "R2 paired must be empty");
        assertTrue(readSequences(r1u).isEmpty(), "R1 unpaired must be empty, tech never goes to unpaired");
        assertTrue(readSequences(r2u).isEmpty(), "R2 unpaired must be empty");
    }

    @Test
    public void testTechRead1_techReadNeverInUnpairedOutput() throws Exception {
        // Deliberately verify the critical invariant: tech read NOT in unpaired even when bio drops.
        File r1in = writeFastq("r1.fastq",
                fqRecord("p1/1", "ACGTACGTACGTACGT"),   // 16 bp tech, long enough to survive if trimmed
                fqRecord("p2/1", "ACGTACGTACGTACGT")); // second pair
        File r2in = writeFastq("r2.fastq",
                fqRecord("p1/2", "ACG"),  // 3 bp bio → dropped
                fqRecord("p2/2", "ACG")); // 3 bp bio → dropped
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(10));

        assertTrue(readSequences(r1u).isEmpty(),
                "Technical read MUST NEVER appear in the unpaired output");
        assertTrue(readSequences(r2u).isEmpty());
        assertTrue(readSequences(r1p).isEmpty());
        assertTrue(readSequences(r2p).isEmpty());
    }

    // -----------------------------------------------------------------------
    // technicalRead=1, the tech read itself is shorter than MinLen threshold

    @Test
    public void testTechRead1_shortTechReadIsPreserved() throws Exception {
        // R1 is only 4bp, MinLen(10) would drop it if applied. Since R1 is the tech
        // read, the trimmer must NOT touch it, and the pair must survive.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGT"));              // 4 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGTACGTACGTACGT")); // 16 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(10));

        assertEquals(List.of("ACGT"),              readSequences(r1p),
                "4bp tech read must survive (MinLen not applied to it)");
        assertEquals(List.of("ACGTACGTACGTACGT"), readSequences(r2p));
        assertTrue(readSequences(r1u).isEmpty());
        assertTrue(readSequences(r2u).isEmpty());
    }

    // -----------------------------------------------------------------------
    // technicalRead=2 (R2=tech, R1=bio)

    @Test
    public void testTechRead2_bioSurvives_bothInPairedOutput() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGTACGTACGT")); // 16 bp bio
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGTACGT"));          // 8 bp tech
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 2, new MinLenTrimmer(5));

        assertEquals(List.of("ACGTACGTACGTACGT"), readSequences(r1p));
        assertEquals(List.of("ACGTACGT"),          readSequences(r2p));
        assertTrue(readSequences(r1u).isEmpty());
        assertTrue(readSequences(r2u).isEmpty());
    }

    @Test
    public void testTechRead2_bioDropped_bothDiscarded() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGT"));              // 4 bp bio
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGTACGTACGTACGT")); // 16 bp tech
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 2, new MinLenTrimmer(10));

        assertTrue(readSequences(r1p).isEmpty());
        assertTrue(readSequences(r2p).isEmpty());
        assertTrue(readSequences(r1u).isEmpty(),
                "Bio read must not go to unpaired when bio is dropped in techRead=2 mode");
        assertTrue(readSequences(r2u).isEmpty(),
                "Tech read must never go to unpaired");
    }

    @Test
    public void testTechRead2_techReadUnchangedByHeadCrop() throws Exception {
        // R2 (tech) has N's, HeadCrop(4) would strip them if applied; must not be applied.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "NNNNTTTTTTTT"));  // 12 bp bio
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "NNNNACGT"));       // 8 bp tech
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 2, new HeadCropTrimmer("4"));

        assertEquals(List.of("TTTTTTTT"), readSequences(r1p), "Bio read must be head-cropped");
        assertEquals(List.of("NNNNACGT"), readSequences(r2p), "Tech read must be unchanged");
    }

    // -----------------------------------------------------------------------
    // Multiple pairs, mixed bio outcomes

    @Test
    public void testTechRead1_multiplePairs_mixedOutcomes() throws Exception {
        // Pair 1: bio=16bp → survives  (MinLen 10)
        // Pair 2: bio= 4bp → dropped
        // Pair 3: bio=20bp → survives
        File r1in = writeFastq("r1.fastq",
                fqRecord("p1/1", "ACGTACGTACGT"),          // 12 bp tech
                fqRecord("p2/1", "ACGTACGTACGT"),          // 12 bp tech
                fqRecord("p3/1", "ACGTACGTACGT"));         // 12 bp tech
        File r2in = writeFastq("r2.fastq",
                fqRecord("p1/2", "ACGTACGTACGTACGT"),      // 16 bp bio → survives
                fqRecord("p2/2", "ACGT"),                  //  4 bp bio → dropped
                fqRecord("p3/2", "ACGTACGTACGTACGTACGT")); // 20 bp bio → survives
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(10));

        List<String> r1Paired = readSequences(r1p);
        List<String> r2Paired = readSequences(r2p);

        assertEquals(2, r1Paired.size(),  "2 of 3 pairs should survive (tech side)");
        assertEquals(2, r2Paired.size(),  "2 of 3 pairs should survive (bio side)");
        assertEquals("ACGTACGTACGT",           r1Paired.get(0));
        assertEquals("ACGTACGTACGT",           r1Paired.get(1));
        assertEquals("ACGTACGTACGTACGT",       r2Paired.get(0));
        assertEquals("ACGTACGTACGTACGTACGT",   r2Paired.get(1));
        assertTrue(readSequences(r1u).isEmpty(), "Tech read must never go to unpaired");
        assertTrue(readSequences(r2u).isEmpty(), "Dropped bio read must not go to unpaired");
    }

    // -----------------------------------------------------------------------
    // Normal PE (technicalRead=0) must be unaffected

    @Test
    public void testNormalPE_techRead0_unpairedOutputStillWorks() throws Exception {
        // In normal PE: R1 short (dropped), R2 long (survives) → R2 goes to R2 unpaired.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGT"));              // 4 bp
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "ACGTACGTACGTACGT")); // 16 bp
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 0, new MinLenTrimmer(10));

        assertTrue(readSequences(r1p).isEmpty(),  "No paired output expected");
        assertTrue(readSequences(r1u).isEmpty(),  "R1 too short, not in unpaired either");
        assertTrue(readSequences(r2p).isEmpty(),  "No paired output expected");
        assertEquals(List.of("ACGTACGTACGTACGT"), readSequences(r2u),
                "Long R2 must be in R2 unpaired in normal PE mode");
    }

    // -----------------------------------------------------------------------
    // Multiple trimmers compose correctly on the bio read

    @Test
    public void testTechRead1_multipleTrimmersCompose() throws Exception {
        // HeadCrop(4) then MinLen(6) on bio read.
        // R2 = "NNNNTTTTTTTT" (12bp) → after HeadCrop(4) = "TTTTTTTT" (8bp) → MinLen(6): 8>=6, survives.
        // R1 (tech): must be unchanged.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "NNNNACGT"));         // 8 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "NNNNTTTTTTTT"));     // 12 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1,
                new HeadCropTrimmer("4"), new MinLenTrimmer(6));

        assertEquals(List.of("NNNNACGT"),   readSequences(r1p), "Tech read unchanged");
        assertEquals(List.of("TTTTTTTT"),   readSequences(r2p), "Bio read head-cropped to 8bp");
        assertTrue(readSequences(r1u).isEmpty());
        assertTrue(readSequences(r2u).isEmpty());
    }

    @Test
    public void testTechRead1_multipleTrimmers_bioDroppedBySecond() throws Exception {
        // HeadCrop(8) leaves bio with 4bp → MinLen(6) drops it → both discarded.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGTACGT"));     // 12 bp tech
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "NNNNNNNNTTTT"));     // 12 bp bio
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u, 1,
                new HeadCropTrimmer("8"), new MinLenTrimmer(6)); // 12-8=4 < 6 → dropped

        assertTrue(readSequences(r1p).isEmpty());
        assertTrue(readSequences(r2p).isEmpty());
        assertTrue(readSequences(r1u).isEmpty(), "Tech read must not appear in unpaired");
        assertTrue(readSequences(r2u).isEmpty());
    }

    // -----------------------------------------------------------------------
    // Empty input

    @Test
    public void testTechRead1_emptyInput_noException() throws Exception {
        File r1in = writeFastq("r1.fastq");  // empty
        File r2in = writeFastq("r2.fastq");  // empty
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        // Must complete without exception, all outputs empty.
        run(r1in, r2in, r1p, r1u, r2p, r2u, 1, new MinLenTrimmer(10));

        assertTrue(readSequences(r1p).isEmpty());
        assertTrue(readSequences(r1u).isEmpty());
        assertTrue(readSequences(r2p).isEmpty());
        assertTrue(readSequences(r2u).isEmpty());
    }

    // -----------------------------------------------------------------------
    // Flag parsing, invalid value via run()

    @Test
    public void testFlagParsing_invalidTechReadValue_returnsFalse() throws Exception {
        // -technicalread 3 is out of range: must be 1 or 2.
        // Provide dummy filenames so the arg count check passes, run() returns false
        // before trying to open any files.
        String[] args = {
            "-technicalread", "3",
            "-phred33",
            "r1.fastq", "r2.fastq",
            "r1p.fastq", "r1u.fastq", "r2p.fastq", "r2u.fastq",
            "MINLEN:10"
        };
        boolean result = TrimmomaticPE.run(args);
        assertFalse(result, "run() must return false for -technicalread 3");
    }
}
