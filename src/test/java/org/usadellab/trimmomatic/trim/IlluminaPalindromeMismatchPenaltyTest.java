package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

/**
 * Palindrome mode must charge Q/10 per mismatch (as simple mode and the README say), not
 * the integer part of Q/10: a mismatch at Q19 costs 1.9, not 1.
 */
public class IlluminaPalindromeMismatchPenaltyTest {
	// TruSeq3-PE prefixes
	private static final String P1 = "TACACTCTTTCCCTACACGACGCTCTTCCGATCT";
	private static final String P2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

	// 2x50 pair whose 42-base insert reads through into the adapters; read 1 carries two
	// sequencing errors (positions 3 and 8) at Q19 ('4'), everything else is Q40 ('I').
	private static final String R1 = "CAGCTTTTTATATTATGCAGAAAATCTACTTCGCCTGATACGAGATCGGA";
	private static final String Q1 = "III4IIII4IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
	private static final String R2 = "CGTATCAGGCGAAGTAGATTTTCTGCATAATATGAAAATCTGAGATCGGA";
	private static final String Q2 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";

	private IlluminaClippingTrimmer makeTrimmer() {
		// seed mismatches 2, palindrome threshold 30, simple threshold 10, min adapter length 8
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(mock(Logger.class), 2, 30, 10, 8, false);
		trimmer.addPrefixPair(P1, P2);
		return trimmer;
	}

	@Test
	public void testTwoMismatchesAtQ19StayBelowThreshold() {
		// 58 aligned bases, 56 matches: 56 * 0.60206 - 2 * 1.9 = 29.9 < 30 -> no palindrome clip.
		// (With the integer penalty the score is 56 * 0.60206 - 2 * 1 = 31.7 and the pair is clipped.)
		FastqRecord[] result = makeTrimmer().processRecords(new FastqRecord[] {
				new FastqRecord("r1", R1, "", Q1, 33), new FastqRecord("r2", R2, "", Q2, 33) });

		assertNotNull(result[0]);
		assertNotNull(result[1]);
		assertEquals(50, result[0].getLength());
		assertEquals(50, result[1].getLength());
	}

	@Test
	public void testOneMismatchAtQ19ClipsThePair() {
		// Revert the second error: 57 matches, 57 * 0.60206 - 1.9 = 32.4 >= 30 -> forward read
		// clipped to the 42-base insert, reverse read dropped (keepBothReads=false).
		String r1 = R1.substring(0, 8) + "C" + R1.substring(9);
		String q1 = Q1.substring(0, 8) + "I" + Q1.substring(9);
		FastqRecord[] result = makeTrimmer().processRecords(new FastqRecord[] {
				new FastqRecord("r1", r1, "", q1, 33), new FastqRecord("r2", R2, "", Q2, 33) });

		assertNotNull(result[0]);
		assertEquals(42, result[0].getLength());
		assertNull(result[1]);
	}
}
