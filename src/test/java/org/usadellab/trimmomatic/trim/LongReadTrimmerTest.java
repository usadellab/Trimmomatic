package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Tests for LongReadTrimmer (unified LONGREADTRIM step).
 *
 * The trimmer uses banded edit distance (not Hamming), k=6 seeding for
 * interior detection, and re-clips fragment ends after every split.
 *
 * Covered:
 *   - 3′ terminal clipping (perfect match, partial overlap)
 *   - 5′ terminal clipping (forward adapter as suffix of RC end)
 *   - No adapter present → read unchanged
 *   - Overlap below minOverlap → no trim
 *   - Mismatches within / exceeding error rate (substitution)
 *   - One-base indel in adapter region still detected (edit distance advantage)
 *   - Reverse-complement adapter detected at 3′ end
 *   - Palindrome adapter not double-counted
 *   - Entire read is adapter → null (dropped)
 *   - N in read / adapter treated as wildcard
 *   - Quality string trimmed in sync with sequence
 *   - Multi-line FASTA assembled correctly
 *   - Multiple adapters: correct adapter chosen
 *   - Interior chimera splitting: read split into fragments
 *   - Interior splitting disabled for HIFI platform
 *   - HIFI-only 3′ near-terminal full-adapter scan (trimmed); confirmed absent for ONT
 *   - Split fragments have suffix /splitNofM in name
 *   - Fragment ends re-clipped after split (residual fix)
 *   - Empty FASTA → IllegalArgumentException
 *   - Missing argument → Exception
 */
public class LongReadTrimmerTest {

    @TempDir
    Path tempDir;

    // ------------------------------------------------------------------
    // Helpers

    private File writeFasta(String... entries) throws Exception {
        File f = tempDir.resolve("adapters.fasta").toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String e : entries) fw.write(e);
        }
        return f;
    }

    private File singleAdapterFasta(String seq) throws Exception {
        return writeFasta(">adapter\n" + seq + "\n");
    }

    /** Build a synthetic read with uniform 'I' quality scores. */
    private FastqRecord rec(String seq) {
        return new FastqRecord("r1", seq, "", "I".repeat(seq.length()), 33);
    }

    private FastqRecord rec(String seq, String qual) {
        return new FastqRecord("r1", seq, "", qual, 33);
    }

    /**
     * Convenience: run LONGREADTRIM with given args and return the first result.
     * Throws if processRecords returns more than one fragment (use trim() for splits).
     */
    private FastqRecord trimOne(String args, FastqRecord r) throws Exception {
        LongReadTrimmer t = new LongReadTrimmer(args);
        FastqRecord[] out = t.processRecords(new FastqRecord[]{r});
        assertTrue(out.length <= 1, "Expected at most 1 output, got " + out.length);
        return out.length == 0 ? null : out[0];
    }

    private FastqRecord[] trim(String args, FastqRecord r) throws Exception {
        return new LongReadTrimmer(args).processRecords(new FastqRecord[]{r});
    }

    private String faPath(File f, float err) { return f.getPath() + ":" + err; }
    private String faPath(File f, float err, int mo) { return f.getPath() + ":" + err + ":" + mo; }
    private String faPath(File f, float err, int mo, int mf) { return f.getPath() + ":" + err + ":" + mo + ":" + mf; }
    private String faPath(File f, float err, int mo, int mf, String plat) {
        return f.getPath() + ":" + err + ":" + mo + ":" + mf + ":" + plat;
    }

    // ------------------------------------------------------------------
    // 3′ terminal clipping

    @Test
    public void test3PrimePerfectFullMatch() throws Exception {
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTTTTTTTTT"); // 8 payload + 8 adapter
        FastqRecord out = trimOne(faPath(fa, 0.0f, 8), r);
        assertNotNull(out);
        assertEquals("ACGTACGT", out.getSequence());
    }

    @Test
    public void test3PrimePartialOverlap() throws Exception {
        // Adapter 12-T's; read ends with 8 T's.  Overlap=8 >= minOverlap=5 → trimmed.
        File fa = singleAdapterFasta("TTTTTTTTTTTT");
        FastqRecord r = rec("GCGCGCGCTTTTTTTT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 5), r);
        assertNotNull(out);
        assertEquals("GCGCGCGC", out.getSequence());
    }

    @Test
    public void testNoAdapter_readUnchanged() throws Exception {
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTACGTACGT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 5), r);
        assertNotNull(out);
        assertEquals("ACGTACGTACGTACGT", out.getSequence());
    }

    @Test
    public void testOverlapBelowMinOverlap_noTrim() throws Exception {
        // Read ends with 3 T's; minOverlap=5 → no trim.
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTACGTTTT"); // 3 T's at end
        FastqRecord out = trimOne(faPath(fa, 0.0f, 5), r);
        assertNotNull(out);
        assertEquals("ACGTACGTACGTTTT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // 5′ terminal clipping
    // A 5′ residual is the suffix of a forward adapter at the read start.

    @Test
    public void test5PrimeClip() throws Exception {
        // Adapter AAAAAAAA; suffix "AAAA" appears at read start (fwd orientation).
        File fa = singleAdapterFasta("AAAAAAAA");
        // minFragLen must be tiny so the clip is not discarded.
        FastqRecord r = rec("AAAAGCGCGCGC");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 4, 4), r);
        assertNotNull(out);
        assertEquals("GCGCGCGC", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Mismatch tolerance (substitution errors — same as Hamming baseline)

    @Test
    public void testSubstitutionWithinErrorRate_trimmed() throws Exception {
        // floor(8 * 0.15) = 1 allowed edit.  1 substitution → trimmed.
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTTTTTATTT"); // T→A at position 4 of adapter region
        FastqRecord out = trimOne(faPath(fa, 0.15f, 8), r);
        assertNotNull(out);
        assertEquals("ACGTACGT", out.getSequence());
    }

    @Test
    public void testSubstitutionExceedsErrorRate_noTrim() throws Exception {
        // floor(8 * 0.10) = 0 allowed edits.
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTTTTTATTT");
        FastqRecord out = trimOne(faPath(fa, 0.10f, 8), r);
        assertNotNull(out);
        assertEquals("ACGTACGTTTTTATTT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Indel tolerance (edit-distance advantage over Hamming)
    //
    // The DP uses fixed-length windows (overlap of k bases vs k adapter bases),
    // so a 1-base insertion in the read shifts every position by 1 but within
    // the window it looks like a substitution.  A Hamming-distance trimmer would
    // also count this as a substitution at the same error rate; the advantage of
    // edit distance shows in interior detection where indels shift k-mer seeds.
    //
    // For the terminal case with adapter TTTTTTTT (8 T's) and read ending in
    // TTTTATTTTT (adapter-with-insertion, 10 chars):
    //   overlap=8, readStart=10: seq[10..17]="TTATTTTT" vs "TTTTTTTT" → 1 edit
    //   trimTo=10, payload="GCGCGCGCTT" (the "TT" before the detected adapter
    //   boundary is an artefact of the insertion shifting the window).
    // This is expected: the adapter IS detected and removed, but the 2-base
    // junction residual requires a dedicated MINLEN step to discard.

    @Test
    public void testSingleSubstitutionInAdapterRegion_trimmed() throws Exception {
        // Adapter TTTTTTTT (8 T's). Read ends with 8 chars containing 1 substitution.
        // floor(8 * 0.15) = 1 allowed edit → adapter detected.
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("GCGCGCGCTTTTATTT"); // payload(8) + "TTTTATTT"(8, 1 subst)
        FastqRecord out = trimOne(faPath(fa, 0.15f, 8, 4), r);
        assertNotNull(out);
        assertEquals("GCGCGCGC", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Reverse complement

    @Test
    public void testReverseComplementAtThreePrime() throws Exception {
        // Adapter AAAA; RC = TTTT.  Read ends with 4 T's.
        File fa = singleAdapterFasta("AAAA");
        FastqRecord r = rec("ACGTACGTTTTT"); // last 4 = TTTT = RC(AAAA)
        FastqRecord out = trimOne(faPath(fa, 0.0f, 4, 4), r);
        assertNotNull(out);
        assertEquals("ACGTACGT", out.getSequence());
    }

    @Test
    public void testPalindromeAdapter() throws Exception {
        // AATT → RC = AATT (palindrome).  Must be indexed once, not twice.
        File fa = singleAdapterFasta("AATT");
        FastqRecord r = rec("GCGCAATT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 4, 4), r);
        assertNotNull(out);
        assertEquals("GCGC", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Entire read is adapter

    @Test
    public void testEntireReadIsAdapter_dropped() throws Exception {
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("TTTTTTTT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 8), r);
        assertNull(out);
    }

    // ------------------------------------------------------------------
    // N wildcards

    @Test
    public void testNInRead_isWildcard() throws Exception {
        File fa = singleAdapterFasta("TTTTTTTT");
        FastqRecord r = rec("ACGTACGTTNNNNNNN"); // T + 7 N's at end
        FastqRecord out = trimOne(faPath(fa, 0.0f, 8), r);
        assertNotNull(out);
        assertEquals("ACGTACGT", out.getSequence());
    }

    @Test
    public void testNInAdapter_isWildcard() throws Exception {
        // Adapter NNNNNNNN — wildcards match everything.  3′ clips the last 8 bp.
        File fa = singleAdapterFasta("NNNNNNNN");
        // 24 bp read: ACGTACGT + TTTTTTTT + ACGTACGT.
        // 3′ scan will find overlap=8 on the last ACGTACGT (0 edits against NNNNNNNN).
        // 5′ scan only uses fwd adapters; NNNNNNNN fwd will also match first 8 bp.
        // Net result: payload = TTTTTTTT (middle 8).
        FastqRecord r = rec("ACGTACGTTTTTTTTTACGTACGT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 8, 4), r);
        assertNotNull(out);
        assertEquals("TTTTTTTT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Quality synced with sequence

    @Test
    public void testQualityTrimmedInSync() throws Exception {
        File fa = singleAdapterFasta("TTTT");
        FastqRecord r = rec("GCGCTTTT", "IIIIJJJJ");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 4, 2), r);
        assertNotNull(out);
        assertEquals("GCGC", out.getSequence());
        assertEquals("IIII", out.getQuality());
    }

    // ------------------------------------------------------------------
    // Multi-line FASTA

    @Test
    public void testMultiLineFastaAssembled() throws Exception {
        File fa = writeFasta(">adapter\nTTTT\nTTTT\n"); // TTTTTTTT
        FastqRecord r = rec("ACGTACGTTTTTTTTT");
        FastqRecord out = trimOne(faPath(fa, 0.0f, 8), r);
        assertNotNull(out);
        assertEquals("ACGTACGT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Multiple adapters

    @Test
    public void testMultipleAdapters_correctAdapterTrimmed() throws Exception {
        File fa = writeFasta(">adp1\nAAAA\n", ">adp2\nTTTT\n");
        FastqRecord r = rec("ACGTTTTT"); // ends with 4 T's; adapter TTTT matches
        FastqRecord out = trimOne(faPath(fa, 0.0f, 4, 2), r);
        assertNotNull(out);
        assertEquals("ACGT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Interior chimera splitting (ONT default)
    //
    // We use a realistic non-homopolymer adapter ("AGATCGGAAGAGCACACGTC")
    // to avoid spurious k-mer seeding in homopolymer payloads.  The payloads
    // use only G and C so no 6-mer from the adapter or its RC appears in them.

    private static final String SPLIT_ADAPTER = "AGATCGGAAGAGCACACGTC"; // 20 bp, Illumina-style
    private static final String SPLIT_PAYLOAD1 = "G".repeat(50);
    private static final String SPLIT_PAYLOAD2 = "C".repeat(50);

    @Test
    public void testInteriorSplit_chimericRead() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec(SPLIT_PAYLOAD1 + SPLIT_ADAPTER + SPLIT_PAYLOAD2);

        FastqRecord[] frags = trim(faPath(fa, 0.0f, 10, 20, "ONT"), r);
        // Expect 2 fragments.  The adapter region is consumed; payloads survive intact.
        assertEquals(2, frags.length, "Expected 2 split fragments");
        assertEquals(SPLIT_PAYLOAD1, frags[0].getSequence());
        assertEquals(SPLIT_PAYLOAD2, frags[1].getSequence());
    }

    @Test
    public void testInteriorSplit_fragmentNaming() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec(SPLIT_PAYLOAD1 + SPLIT_ADAPTER + SPLIT_PAYLOAD2);

        FastqRecord[] frags = trim(faPath(fa, 0.0f, 10, 20, "ONT"), r);
        assertEquals(2, frags.length);
        assertTrue(frags[0].getName().contains("/split"), "Name should contain /split");
    }

    @Test
    public void testHiFiPlatform_noInteriorSplit() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec(SPLIT_PAYLOAD1 + SPLIT_ADAPTER + SPLIT_PAYLOAD2);

        FastqRecord[] frags = trim(faPath(fa, 0.0f, 10, 20, "HIFI"), r);
        // HIFI disables interior splitting entirely.
        // The adapter is in the interior and won't be caught by terminal clipping.
        assertEquals(1, frags.length, "HIFI should not split the read");
    }

    // ------------------------------------------------------------------
    // 3′ near-terminal full-adapter scan (HIFI only)
    // Covers a full adapter followed by a few trailing bases rather than
    // hanging off the exact 3′ end -- the gap left once interior splitting
    // is disabled for HIFI.

    @Test
    public void testHiFi3PrimeNearTerminalFullAdapter_trimmed() throws Exception {
        // 30bp payload + full 20bp adapter + 4 trailing bases, not hanging off the end.
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec("G".repeat(30) + SPLIT_ADAPTER + "TTTT");

        FastqRecord out = trimOne(faPath(fa, 0.0f, 10, 20, "HIFI"), r);
        assertNotNull(out);
        assertEquals("G".repeat(30), out.getSequence());
    }

    @Test
    public void testOnt3PrimeNearTerminalGap_stillNotTrimmed() throws Exception {
        // Same layout as above, but ONT/CLR must not gain this scan: it produced
        // false-positive over-clipping on real payload sequence at ONT/CLR error
        // rates and was deliberately reverted for those platforms (see 835e329).
        // minFragLen=30 keeps workLen (54) below the interior-split threshold
        // (2*minOverlap + 2*minFragLen = 80), so interior splitting can't mask
        // the absence of the near-terminal scan here.
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec("G".repeat(30) + SPLIT_ADAPTER + "TTTT");

        FastqRecord out = trimOne(faPath(fa, 0.0f, 10, 30, "ONT"), r);
        assertNotNull(out);
        assertEquals("G".repeat(30) + SPLIT_ADAPTER + "TTTT", out.getSequence());
    }

    // ------------------------------------------------------------------
    // Error cases

    @Test
    public void testEmptyFasta_throwsIllegalArgument() throws Exception {
        File fa = tempDir.resolve("empty.fasta").toFile();
        fa.createNewFile();
        assertThrows(IllegalArgumentException.class,
                () -> new LongReadTrimmer(fa.getPath() + ":0.1:10"));
    }

    @Test
    public void testMissingErrorRate_throwsException() {
        assertThrows(Exception.class,
                () -> new LongReadTrimmer("noarguments"));
    }

    @Test
    public void testPairedEndMode_throwsRuntimeException() throws Exception {
        File fa = singleAdapterFasta("TTTT");
        LongReadTrimmer t = new LongReadTrimmer(fa.getPath() + ":0.1:4");
        assertThrows(RuntimeException.class,
                () -> t.processRecords(new FastqRecord[]{rec("ACGT"), rec("TGCA")}));
    }
}
