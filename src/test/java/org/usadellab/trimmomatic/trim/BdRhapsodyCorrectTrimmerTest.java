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

public class BdRhapsodyCorrectTrimmerTest {

    // Arbitrary but distinct 9bp codes -- the trimmer never validates L1/L2
    // content, only their length, so filler bases can be anything.
    private static final String CLS1 = "AAACCCGGT";
    private static final String CLS2 = "TTTGGGCCA";
    private static final String CLS3 = "GGGTTTAAC";
    private static final String UMI = "ACGTACGT";
    private static final String L1_FILLER = "NNNNNNNNNNNN";    // 12bp, nominal (V1)
    private static final String L2_FILLER = "NNNNNNNNNNNNN";   // 13bp, nominal (V1)

    // ENHANCEDV2 mode's linkers are much shorter than V1's -- 4bp each.
    private static final String L1_FILLER_ENHANCEDV2 = "NNNN";
    private static final String L2_FILLER_ENHANCEDV2 = "NNNN";

    @TempDir
    Path tempDir;

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    private String nominalEnhancedV2Read() {
        return CLS1 + L1_FILLER_ENHANCEDV2 + CLS2 + L2_FILLER_ENHANCEDV2 + CLS3 + UMI;
    }

    private String nominalRead() {
        return CLS1 + L1_FILLER + CLS2 + L2_FILLER + CLS3 + UMI;
    }

    private File writeWhitelistDir(String cls1Entry, String cls2Entry, String cls3Entry) throws IOException {
        File dir = tempDir.resolve("wl_" + System.identityHashCode(new Object())).toFile();
        dir.mkdirs();
        writeFile(new File(dir, "CLS1.txt"), cls1Entry);
        writeFile(new File(dir, "CLS2.txt"), cls2Entry);
        writeFile(new File(dir, "CLS3.txt"), cls3Entry);
        return dir;
    }

    private void writeFile(File f, String content) throws IOException {
        try (FileWriter fw = new FileWriter(f)) {
            fw.write(content + "\n");
        }
    }

    // ------------------------------------------------------------------
    // Nominal (no offset), exact match

    @Test
    public void testNominalRead_exactMatch() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        FastqRecord result = trimmer.processRecord(makeRecord("r1", nominalRead()));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
        assertEquals("", result.getSequence());
    }

    @Test
    public void testPayloadAfterUmi_survives() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        FastqRecord result = trimmer.processRecord(makeRecord("r1", nominalRead() + "TTTTT"));

        assertNotNull(result);
        assertEquals("TTTTT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Offset handling via CLS2

    @Test
    public void testOneBaseDeletion_offsetMinusOne_stillLocatesCls3AndUmi() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String shortL1 = L1_FILLER.substring(1); // 11bp instead of 12 -- 1bp deletion
        // A real read stays at its fixed cycle length regardless of an internal
        // deletion -- the missing base just means one extra trailing base (here,
        // arbitrary carryover) gets read instead of being cut off. Pad back out
        // so the read still clears the 60bp minimum-length check.
        String seq = CLS1 + shortL1 + CLS2 + L2_FILLER + CLS3 + UMI + "A";
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testOneBaseInsertion_offsetPlusOne_stillLocatesCls3AndUmi() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String longL1 = L1_FILLER + "N"; // 13bp instead of 12 -- 1bp insertion
        String seq = CLS1 + longL1 + CLS2 + L2_FILLER + CLS3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testTwoBaseDeletion_offsetMinusTwo_stillLocatesCls3AndUmi() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String shortL1 = L1_FILLER.substring(2); // 10bp instead of 12 -- 2bp deletion
        String seq = CLS1 + shortL1 + CLS2 + L2_FILLER + CLS3 + UMI + "AA"; // pad back to 60bp, see comment above
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testTwoBaseInsertion_offsetPlusTwo_stillLocatesCls3AndUmi() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String longL1 = L1_FILLER + "NN"; // 14bp instead of 12 -- 2bp insertion
        String seq = CLS1 + longL1 + CLS2 + L2_FILLER + CLS3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    // ------------------------------------------------------------------
    // Hamming-1 correction (substitution noise, no offset)

    @Test
    public void testCls1OneMismatch_corrects() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls1 = "TAACCCGGT"; // position 0: A -> T
        String seq = rawCls1 + L1_FILLER + CLS2 + L2_FILLER + CLS3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName(),
                "corrected CLS1 in the name must be the WHITELIST entry, not the raw sequence");
    }

    @Test
    public void testCls2OneMismatch_noOffsetNeeded_correctsAtNominal() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls2 = "ATTGGGCCA"; // position 0: T -> A
        String seq = CLS1 + L1_FILLER + rawCls2 + L2_FILLER + CLS3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testCls3OneMismatch_corrects() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls3 = "AGGTTTAAC"; // position 0: G -> A
        String seq = CLS1 + L1_FILLER + CLS2 + L2_FILLER + rawCls3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    // ------------------------------------------------------------------
    // Drop cases

    @Test
    public void testCls1TwoMismatches_beyondTolerance_drops() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls1 = "TTACCCGGT"; // positions 0 and 1 both wrong
        String seq = rawCls1 + L1_FILLER + CLS2 + L2_FILLER + CLS3 + UMI;
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));
    }

    @Test
    public void testCls2Unresolvable_dropsEvenThoughCls1AndCls3Match() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls2 = "AATGGGCCA"; // two mismatches, no offset explains it either
        String seq = CLS1 + L1_FILLER + rawCls2 + L2_FILLER + CLS3 + UMI;
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));
    }

    @Test
    public void testCls3TwoMismatches_beyondTolerance_drops() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String rawCls3 = "AAGTTTAAC"; // positions 0 and 1 both wrong
        String seq = CLS1 + L1_FILLER + CLS2 + L2_FILLER + rawCls3 + UMI;
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));
    }

    @Test
    public void testShorterThan60bp_drops() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        String seq = nominalRead().substring(0, 59); // one short of the full construct
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));
    }

    // ------------------------------------------------------------------
    // Argument validation

    @Test
    public void testMaxMismatchZero_onlyExactMatchWorks() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":0");

        String rawCls1 = "TAACCCGGT"; // 1 mismatch -- should NOT correct at maxMismatch=0
        String seq = rawCls1 + L1_FILLER + CLS2 + L2_FILLER + CLS3 + UMI;
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));

        assertNotNull(trimmer.processRecord(makeRecord("r2", nominalRead())));
    }

    @Test
    public void testMismatchAboveSupportedCap_throws() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        assertThrows(IllegalArgumentException.class,
                () -> new BdRhapsodyCorrectTrimmer(wl.getPath() + ":2"));
    }

    @Test
    public void testMissingRequiredArgs_throws() throws Exception {
        assertThrows(IllegalArgumentException.class,
                () -> new BdRhapsodyCorrectTrimmer("onlyOneToken"));
    }

    @Test
    public void testCustomSeparator() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1:__");

        FastqRecord result = trimmer.processRecord(makeRecord("r1", nominalRead()));

        assertNotNull(result);
        assertEquals("r1__" + CLS1 + CLS2 + CLS3 + "__" + UMI, result.getName());
    }

    @Test
    public void testMissingWhitelistFile_throws() throws Exception {
        File dir = tempDir.resolve("incomplete_wl").toFile();
        dir.mkdirs();
        writeFile(new File(dir, "CLS1.txt"), CLS1);
        // CLS2.txt / CLS3.txt intentionally missing
        assertThrows(IOException.class, () -> new BdRhapsodyCorrectTrimmer(dir.getPath() + ":1"));
    }

    // ------------------------------------------------------------------
    // Bead version selection (V1 default / explicit / ENHANCEDV2)

    @Test
    public void testExplicitV1Token_behavesSameAsDefault() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("V1:" + wl.getPath() + ":1");

        FastqRecord result = trimmer.processRecord(makeRecord("r1", nominalRead()));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testEnhancedV2_zeroInset_exactMatch() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        FastqRecord result = trimmer.processRecord(makeRecord("r1", nominalEnhancedV2Read()));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
        assertEquals("", result.getSequence());
    }

    @Test
    public void testEnhancedV2_insetOfOne_locatesCorrectly() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        String seq = "A" + nominalEnhancedV2Read(); // 1bp prefix inset
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testEnhancedV2_insetOfTwo_locatesCorrectly() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        String seq = "GT" + nominalEnhancedV2Read(); // 2bp prefix inset
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testEnhancedV2_insetOfThree_locatesCorrectly() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        String seq = "TCA" + nominalEnhancedV2Read(); // 3bp prefix inset
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testEnhancedV2_cls2OneMismatch_correctsAtNominalOffset() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        String rawCls2 = "ATTGGGCCA"; // position 0: T -> A
        String seq = CLS1 + L1_FILLER_ENHANCEDV2 + rawCls2 + L2_FILLER_ENHANCEDV2 + CLS3 + UMI;
        FastqRecord result = trimmer.processRecord(makeRecord("r1", seq));

        assertNotNull(result);
        assertEquals("r1_" + CLS1 + CLS2 + CLS3 + "_" + UMI, result.getName());
    }

    @Test
    public void testEnhancedV2_cls1Unresolvable_drops() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer("ENHANCEDV2:" + wl.getPath() + ":1");

        // No inset candidate (0-3) lines up CLS1 with the whitelist entry, and
        // no offset qualifies as a 1-mismatch neighbour of it either.
        String garbledCls1 = "CCCAAAGGT";
        String seq = garbledCls1 + L1_FILLER_ENHANCEDV2 + CLS2 + L2_FILLER_ENHANCEDV2 + CLS3 + UMI;
        assertNull(trimmer.processRecord(makeRecord("r1", seq)));
    }

    // ------------------------------------------------------------------
    // Symmetric-mode guard

    @Test
    public void testProcessRecords_twoElementArray_throws() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        FastqRecord a = makeRecord("a", nominalRead());
        FastqRecord b = makeRecord("b", nominalRead());
        assertThrows(IllegalStateException.class, () -> trimmer.processRecords(new FastqRecord[] { a, b }));
    }

    @Test
    public void testProcessRecords_singleElementArray_works() throws Exception {
        File wl = writeWhitelistDir(CLS1, CLS2, CLS3);
        BdRhapsodyCorrectTrimmer trimmer = new BdRhapsodyCorrectTrimmer(wl.getPath() + ":1");

        FastqRecord[] in = { makeRecord("r1", nominalRead()) };
        FastqRecord[] result = trimmer.processRecords(in);

        assertEquals(1, result.length);
        assertNotNull(result[0]);
    }
}
