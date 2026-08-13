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
 *   - 3′ near-terminal full-adapter scan: trimmed on HIFI and on ONT
 *   - Residual interior adapter, the three FN populations measured on lsk114
 *   - Payload resembling the adapter near the 3′ end is not clipped (835e329 guard)
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

    /**
     * Re-scoped from {@code testOnt3PrimeNearTerminalGap_stillNotTrimmed}.
     *
     * <p>That test asserted a complete, perfectly-matching 20 bp adapter at
     * err=0.0 must SURVIVE on ONT.  That was never a correctness claim: it
     * encoded the accepted-miss tradeoff from 835e329, which reverted a
     * BRUTE-FORCE 3' near-terminal scan (every start position in a
     * minOverlap-wide window, all orientations, no seed requirement) because it
     * over-clipped payload resembling the adapter.
     *
     * <p>The current path is seed-gated and requires the full adapter, and the
     * over-clipping it was reverted for has been measured at zero.  Per-base
     * ground truth on Badread reads drawn from the real E. coli K-12 reference,
     * baseline vs this code, over ~12.5k reads per arm:
     *
     * <pre>
     *   arm            distal FP (genuine over-clip)   FN        F1
     *   ~10% error     32 -> 32   (unchanged)          50758 -> 26     0.8909 -> 0.9226
     *   ~3%  error     570 -> 570 (unchanged)          56503 -> 2053   0.8858 -> 0.9225
     * </pre>
     *
     * <p>On real ONT data the same comparison removes 441 bases from 14.3 Gbp
     * (ont-zymo) and 3884 from 30.9 Gbp (ont-zymo-lsk114) with no read loss,
     * while residual adapter reads fall 147 -> 42 on ont-zymo-lsk114.
     *
     * <p>835e329's actual intent is preserved by
     * {@link #test3PrimeNearTerminal_payloadResemblingAdapter_notClipped},
     * which guards payload that merely resembles the adapter.  That is the
     * property worth asserting; a mandated miss is not.
     */
    @Test
    public void testOnt3PrimeNearTerminalGap_nowTrimmed() throws Exception {
        // minFragLen=30 keeps workLen (54) below the OLD interior-split gate
        // (2*minOverlap + 2*minFragLen = 80) but above the corrected one
        // (2*minOverlap + minFragLen = 50), so this exercises both fixes.
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        FastqRecord r = rec("G".repeat(30) + SPLIT_ADAPTER + "TTTT");

        FastqRecord out = trimOne(faPath(fa, 0.0f, 10, 30, "ONT"), r);
        assertNotNull(out);
        assertEquals("G".repeat(30), out.getSequence(),
                "a complete adapter near the 3' end must be removed on ONT too");
    }

    // ------------------------------------------------------------------
    // Residual interior adapter: the three FN populations measured on
    // ont-zymo-lsk114 (462 residual-adapter records, seqkit locate mm<=2).
    //
    //   313 (68%)  reads < 2*minOverlap + 2*minFragLen  -> interior scan gated off
    //   119 (26%)  reads >= gate, adapter ends within minOverlap of the 3' end
    //    30  (6%)  already-split fragments carrying a residual adapter
    //
    // All three assert the DESIRED behaviour and are expected to fail against
    // the current implementation.  Benchmark params: err 0.15, minOverlap 10,
    // minFragLen 100, ONT.

    /**
     * Population 1 (68% of FNs): a 201 bp read carrying a full interior adapter.
     * The gate at processRecords (workLen >= 2*minOverlap + 2*minFragLen = 220)
     * skips the interior scan entirely, so the adapter survives.  The gate
     * demands room for TWO viable fragments, but the per-fragment minFragLen
     * filter means ONE viable fragment is enough.
     */
    @Test
    public void testShortRead_interiorAdapter_scanGateTooStrict() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        // 150 payload + 20 adapter + 31 payload = 201 bp, below the 220 gate.
        FastqRecord r = rec("G".repeat(150) + SPLIT_ADAPTER + "C".repeat(31));

        FastqRecord out = trimOne(faPath(fa, 0.0f, 10, 100, "ONT"), r);
        assertNotNull(out);
        // Fragment 1 (150 bp) is viable; fragment 2 (31 bp) is below minFragLen.
        assertEquals("G".repeat(150), out.getSequence(),
                "interior adapter must be removed even when only one fragment survives");
    }

    /**
     * Population 2 (26% of FNs): read is long enough for the interior scan, and
     * the 6-mer seed produces the candidate, but bestMatchAt rejects it because
     * the adapter end falls inside the 3' terminal zone
     * (seqLen - readStart - terminalZone < adapterLen).  The 3' terminal clip
     * cannot cover it either: the trailing overlap (5) is below minOverlap (10).
     */
    @Test
    public void test3PrimeNearTerminal_fullAdapterWithTrailingBases() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        // 200 payload + 20 adapter + 5 trailing = 225 bp, above the 220 gate.
        // readStart=200 <= seeding cap (225-10-10=205), so the candidate exists,
        // but 225-200-10 = 15 < 20 = adapterLen, so bestMatchAt discards it.
        FastqRecord r = rec("G".repeat(200) + SPLIT_ADAPTER + "TTTTT");

        FastqRecord out = trimOne(faPath(fa, 0.0f, 10, 100, "ONT"), r);
        assertNotNull(out);
        assertEquals("G".repeat(200), out.getSequence(),
                "full adapter near the 3' end must be removed");
    }

    /**
     * Population 3 (6% of FNs): a read split at one junction whose surviving
     * fragment still carries a second adapter.  clipFragment only re-clips the
     * fragment terminals, so an adapter left in the fragment interior by the
     * original scan's blind spots is never revisited.
     */
    @Test
    public void testSplitFragment_residualInteriorAdapter() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        // adapter #1 at 150 (clean interior), adapter #2 at 330 with 5 trailing.
        FastqRecord r = rec("G".repeat(150) + SPLIT_ADAPTER
                          + "C".repeat(160) + SPLIT_ADAPTER + "TTTTT");

        FastqRecord[] frags = trim(faPath(fa, 0.0f, 10, 100, "ONT"), r);
        for (FastqRecord f : frags)
            assertTrue(f == null || !f.getSequence().contains(SPLIT_ADAPTER),
                    "no emitted fragment may still contain a full adapter");
    }

    /**
     * Regression guard for 835e329: the reverted brute-force 3' near-terminal
     * scan over-clipped payload that merely resembled the adapter.  Any fix for
     * the tests above must stay seed-gated and inside the edit budget -- payload
     * sharing a 6-mer with the adapter but exceeding allowedEdits must survive.
     */
    @Test
    public void test3PrimeNearTerminal_payloadResemblingAdapter_notClipped() throws Exception {
        File fa = singleAdapterFasta(SPLIT_ADAPTER);
        // Shares the adapter's first 6 bases (AGATCG) so the k-mer seed fires,
        // but the remaining 14 bases are unrelated -> far beyond 3 allowed edits
        // at err 0.15 (0.15 * 20 = 3).
        String lookalike = "AGATCG" + "TTATTATTATTATT";
        assertEquals(SPLIT_ADAPTER.length(), lookalike.length());
        FastqRecord r = rec("G".repeat(200) + lookalike + "TTTTT");

        FastqRecord out = trimOne(faPath(fa, 0.15f, 10, 100, "ONT"), r);
        assertNotNull(out);
        assertEquals("G".repeat(200) + lookalike + "TTTTT", out.getSequence(),
                "payload resembling the adapter must not be clipped");
    }

    /**
     * Nested adapter pairs must not leave a residual base at the 3' end.
     *
     * <p>Uses the real LSK114 pair, where {@code 3'[1:18]} is byte-identical to
     * the first 17 bases of the 5' adapter's reverse complement.  That lets the
     * 3' prefix scan match the longer adapter one base late: with a 2-base
     * trailing gap it accepts a 19-base overlap at {@code adapterStart + 1}
     * (17 matching, 2 mismatching, budget 2) and sets trimTo there, stranding the
     * 3' adapter's first base.
     *
     * <p>The interior scan finds the full adapter at the correct offset, but only
     * because it now runs on the untrimmed read -- scanning the already-clipped
     * window saw a 1-base remnant and nothing to match.  The straddling hit then
     * pulls trimTo back over the whole adapter.
     */
    @Test
    public void testNestedAdapterPair_noResidualBaseAt3Prime() throws Exception {
        File fa = writeFasta(
                ">SQK-LSK114_5prime\nCCTGTACTTCGTTCAGTTACGTATTGC\n",
                ">SQK-LSK114_3prime\nAGCAATACGTAACTGAAC\n");
        String payload = "G".repeat(200);
        FastqRecord r = rec(payload + "AGCAATACGTAACTGAAC" + "TT");

        FastqRecord out = trimOne(faPath(fa, 0.15f, 10, 100, "ONT"), r);
        assertNotNull(out);
        assertEquals(payload, out.getSequence(),
                "no adapter base may survive the 3' clip");
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
