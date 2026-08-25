package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class OrientTrimmerTest {

    /** The scCOLOR-seq SMART primer, the anchor this step was built for. */
    private static final String PRIMER = "AAGCAGTGGTATCAACGCAGAGT";
    private static final String PRIMER_RC = "ACTCTGCGTTGATACCACTGCTT";

    private static Path writeFasta(String... seqs) throws IOException {
        Path p = Files.createTempFile("orient-test", ".fa");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < seqs.length; i++)
            sb.append(">q").append(i).append('\n').append(seqs[i]).append('\n');
        Files.writeString(p, sb.toString());
        p.toFile().deleteOnExit();
        return p;
    }

    private static OrientTrimmer trimmer(String fasta, String rest) throws IOException {
        return new OrientTrimmer(fasta + ":" + rest);
    }

    private FastqRecord makeRecord(String name, String seq) {
        return new FastqRecord(name, seq, "", "I".repeat(seq.length()), 33);
    }

    /** Quality that differs per position, so a reversal is detectable. */
    private FastqRecord makeRecordGradedQuality(String name, String seq) {
        StringBuilder q = new StringBuilder();
        for (int i = 0; i < seq.length(); i++)
            q.append((char) ('!' + (i % 40)));
        return new FastqRecord(name, seq, "", q.toString(), 33);
    }

    // ------------------------------------------------------------------
    // Forward orientation

    @Test
    public void testForwardMatch_passesThroughUnchanged() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        FastqRecord rec = makeRecord("r1", "GGGGGGGGGG" + PRIMER + "ACGTACGTACGT");
        FastqRecord result = t.processRecord(rec);

        assertSame(rec, result, "a forward read should be handed back untouched");
    }

    @Test
    public void testForwardMatchWithEdits_stillForward() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        // Two substitutions inside the primer, within the floor(23 * 0.15) = 3
        // edit budget. Nanopore reads seldom carry the anchor without errors,
        // so exact matching alone would be useless on real data.
        String mutated = PRIMER.substring(0, 5) + "CC" + PRIMER.substring(7);
        assertEquals(PRIMER.length(), mutated.length());

        FastqRecord rec = makeRecord("r1", "TTTTT" + mutated + "ACGTACGT");
        assertSame(rec, t.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Reverse orientation

    @Test
    public void testReverseMatch_readIsReverseComplemented() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        String payload = "ACGTACGTACGT";
        FastqRecord rec = makeRecord("r1", payload + PRIMER_RC + "GGGGGGGGGG");
        FastqRecord result = t.processRecord(rec);

        assertNotNull(result);
        // Reverse-complementing must bring the primer back to forward orientation.
        assertEquals("CCCCCCCCCC" + PRIMER + "ACGTACGTACGT", result.getSequence());
        assertEquals("r1", result.getName(), "orientation must not rename the read");
    }

    @Test
    public void testReverseMatch_qualityIsReversedNotComplemented() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        FastqRecord rec = makeRecordGradedQuality("r1", "ACGT" + PRIMER_RC + "TTTT");
        FastqRecord result = t.processRecord(rec);

        assertNotNull(result);
        String original = rec.getQuality();
        String expected = new StringBuilder(original).reverse().toString();
        assertEquals(expected, result.getQuality());
        assertEquals(result.getSequence().length(), result.getQuality().length());
    }

    // ------------------------------------------------------------------
    // Undecidable reads

    @Test
    public void testNoMatch_droppedByDefault() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        FastqRecord rec = makeRecord("r1", "ACGTACGTACGTACGTACGTACGTACGTACGT");
        assertNull(t.processRecord(rec));
    }

    @Test
    public void testNoMatch_keptWhenOnMissingIsKeep() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15:KEEP");

        FastqRecord rec = makeRecord("r1", "ACGTACGTACGTACGTACGTACGTACGTACGT");
        assertSame(rec, t.processRecord(rec));
    }

    @Test
    public void testBothOrientationsPresent_isTreatedAsUndecidable() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        // A chimera carrying the primer on both strands. The evidence points
        // both ways and there is no basis for choosing one.
        FastqRecord rec = makeRecord("r1", PRIMER + "ACGTACGT" + PRIMER_RC);
        assertNull(t.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Argument handling

    @Test
    public void testTooFewArguments() {
        assertThrows(IllegalArgumentException.class, () -> new OrientTrimmer("only_a_path.fa"));
    }

    @Test
    public void testRejectsErrorRateOutOfRange() throws IOException {
        String fasta = writeFasta(PRIMER).toString();
        assertThrows(IllegalArgumentException.class, () -> trimmer(fasta, "1.5"));
        assertThrows(IllegalArgumentException.class, () -> trimmer(fasta, "-0.1"));
    }

    @Test
    public void testRejectsEmptyFasta() throws IOException {
        Path empty = Files.createTempFile("orient-empty", ".fa");
        empty.toFile().deleteOnExit();
        assertThrows(IllegalArgumentException.class, () -> trimmer(empty.toString(), "0.15"));
    }

    @Test
    public void testRejectsOverlongQuery() throws IOException {
        String tooLong = "A".repeat(65);
        String fasta = writeFasta(tooLong).toString();
        assertThrows(IllegalArgumentException.class, () -> trimmer(fasta, "0.15"));
    }

    @Test
    public void testPathContainingColonIsParsed() throws IOException {
        // Mimics a Windows drive letter: the path must not be split on its colon.
        Path p = writeFasta(PRIMER);
        File f = p.toFile();
        String withColonLikePrefix = f.getAbsolutePath();
        OrientTrimmer t = new OrientTrimmer(withColonLikePrefix + ":0.15");

        FastqRecord rec = makeRecord("r1", "GG" + PRIMER + "ACGT");
        assertNotNull(t.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // Paired-end refusal

    @Test
    public void testPairedEndIsRefused() throws IOException {
        OrientTrimmer t = trimmer(writeFasta(PRIMER).toString(), "0.15");

        FastqRecord[] pair = {
                makeRecord("r1", PRIMER + "ACGT"),
                makeRecord("r1", PRIMER + "ACGT")
        };
        assertThrows(IllegalStateException.class, () -> t.processRecords(pair));
    }

    // ------------------------------------------------------------------
    // Multiple anchors

    @Test
    public void testSecondFastaEntryAlsoOrients() throws IOException {
        String other = "TTTTGGGGCCCCAAAATTTTGG";
        OrientTrimmer t = trimmer(writeFasta(PRIMER, other).toString(), "0.10");

        String otherRc = LongReadTrimmer.reverseComplement(other);
        FastqRecord rec = makeRecord("r1", "ACGT" + otherRc + "ACGT");
        FastqRecord result = t.processRecord(rec);

        assertNotNull(result);
        assertEquals("ACGT" + other + "ACGT", result.getSequence());
    }
}
