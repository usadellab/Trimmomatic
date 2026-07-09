package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	private String makeQuality(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	@Test
	public void testParseValidRecord(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test.fastq");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			FastqRecord rec = parser.next();
			assertEquals("seq1", rec.getName());
			assertEquals(seq, rec.getSequence());
			assertEquals(qual, rec.getQuality());
			
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}
	
	@Test
	public void testParseMultipleRecords(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_multi.fastq");
		String seq1 = makeSequence(150);
		String qual1 = makeQuality('I', 150);
		String seq2 = makeSequence(150);
		String qual2 = makeQuality('J', 150);

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq1 + "\n+\n" + qual1 + "\n");
			writer.write("@seq2\n" + seq2 + "\n+\n" + qual2 + "\n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			assertEquals("seq1", parser.next().getName());
			assertTrue(parser.hasNext());
			assertEquals("seq2", parser.next().getName());
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testInvalidNameLine(@TempDir File tempDir) throws IOException {
		// A line starting with neither '@' nor '>' is truly invalid
		File fastqFile = new File(tempDir, "test_invalid_name.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("Xseq1\nACGT\n+\nIIII\n");
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}

	@Test
	public void testFastaInput(@TempDir File tempDir) throws IOException {
		File fastaFile = new File(tempDir, "test.fasta");
		try (FileWriter writer = new FileWriter(fastaFile)) {
			writer.write(">seq1\nACGTACGT\n>seq2\nTTTT\n");
		}

		FastqParser parser = new FastqParser(33);
		try {
			parser.openFasta(fastaFile);

			assertTrue(parser.hasNext());
			FastqRecord r1 = parser.next();
			assertEquals("seq1", r1.getName());
			assertEquals("ACGTACGT", r1.getSequence());
			assertEquals(8, r1.getQuality().length());

			assertTrue(parser.hasNext());
			FastqRecord r2 = parser.next();
			assertEquals("seq2", r2.getName());
			assertEquals("TTTT", r2.getSequence());

			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testIsFasta(@TempDir File tempDir) throws IOException {
		File fastaFile = new File(tempDir, "test.fasta");
		try (FileWriter writer = new FileWriter(fastaFile)) {
			writer.write(">seq1\nACGT\n");
		}
		assertTrue(FastqParser.isFasta(fastaFile));

		File fastqFile = new File(tempDir, "test.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\nACGT\n+\nIIII\n");
		}
		assertFalse(FastqParser.isFasta(fastqFile));
	}

	@Test
	public void testInvalidCommentLine(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_invalid_comment.fastq");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n-\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}

	@Test
	public void testTruncatedFileMissingQual(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_truncated.fastq");
		String seq = makeSequence(150);
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n");
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}

	@Test
	public void testTruncatedFileMissingSeq(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_truncated_seq.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n"); // Only name line
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}

	@Test
	public void testTruncatedFileMissingComment(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_truncated_comment.fastq");
		String seq = makeSequence(150);
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n"); // Name + Seq, missing +
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}

	@Test
	public void testLengthMismatch(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_mismatch.fastq");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 149); // 1 base shorter
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(33);
		try {
			assertThrows(RuntimeException.class, () -> parser.open(fastqFile));
		} finally {
			parser.close();
		}
	}
	
	@Test
	public void testEmptyFile(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_empty.fastq");
		fastqFile.createNewFile();

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);
		try {
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}
}
