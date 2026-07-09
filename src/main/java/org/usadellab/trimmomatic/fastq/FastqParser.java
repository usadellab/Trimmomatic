package org.usadellab.trimmomatic.fastq;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayDeque;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.util.PositionTrackingInputStream;
import org.usadellab.trimmomatic.util.compression.CompressionFormat;

public class FastqParser implements Closeable {

	private static final int  PREREAD_COUNT     = 10000;
	private static final long PREREAD_MAX_BYTES = 4L * 1024 * 1024; // 4 MB byte cap for phred detection

	private int phredOffset;
	private ArrayDeque<FastqRecord> deque;
	int qualHistogram[];

	private PositionTrackingInputStream posTrackInputStream;
	private BufferedReader reader;
	private FastqRecord current;

	private AtomicBoolean atEOF;
	/** Buffered header line carried over from one FASTA record to the next. */
	private String fastaBufferedHeader = null;

	public FastqParser(int phredOffset) {
		this.phredOffset = phredOffset;
		deque = new ArrayDeque<FastqRecord>(PREREAD_COUNT);

		this.atEOF = new AtomicBoolean();
	}

	public void setPhredOffset(int phredOffset) {
		this.phredOffset = phredOffset;

		if (current != null)
			current.setPhredOffset(phredOffset);
	}

	public void parseOne() throws IOException {
		current = null;

		String line;

		if (fastaBufferedHeader != null) {
			line = fastaBufferedHeader;
			fastaBufferedHeader = null;
		} else {
			line = reader.readLine();
		}

		if (line == null) {
			atEOF.set(true);
			return;
		}

		if (line.charAt(0) == '@') {
			// Standard FASTQ record
			String name = line.substring(1);

			String sequence = reader.readLine();
			if (sequence == null)
				throw new RuntimeException("Missing sequence line from record: " + name);

			String commentLine = reader.readLine();
			if (commentLine == null)
				throw new RuntimeException("Missing comment line from record: " + name);
			if (commentLine.charAt(0) != '+')
				throw new RuntimeException("Invalid FASTQ comment line: " + commentLine);

			String quality = reader.readLine();
			if (quality == null)
				throw new RuntimeException("Missing quality line from record: " + name);

			current = new FastqRecord(name, sequence, commentLine.substring(1), quality, phredOffset);

		} else if (line.charAt(0) == '>') {
			// FASTA record: assemble multi-line sequence, generate dummy Phred+33 quality
			String name = line.substring(1).trim();
			StringBuilder seq = new StringBuilder();
			String seqLine;
			while ((seqLine = reader.readLine()) != null && !seqLine.startsWith(">")) {
				if (!seqLine.startsWith(";"))
					seq.append(seqLine.trim());
			}
			// seqLine is null (EOF) or the next '>' header — buffer it for the next call
			fastaBufferedHeader = seqLine;

			String sequence = seq.toString().toUpperCase();
			String quality = "I".repeat(sequence.length()); // Phred+40, valid Phred+33 dummy
			current = new FastqRecord(name, sequence, "", quality, phredOffset);

		} else {
			throw new RuntimeException("Invalid FASTQ/FASTA name line: " + line);
		}
	}

	public int getProgress() {
		if (atEOF.get())
			return 100;

		return posTrackInputStream.getProgressPercentage();
	}

	private void accumulateHistogram(FastqRecord rec) {
		int quals[] = rec.getQualityAsInteger(false);

		for (int i : quals)
			qualHistogram[i]++;
	}

	public int determinePhredOffset() {
		int phred33Total = 0;
		int phred64Total = 0;

		for (int i = 33; i <= 58; i++)
			phred33Total += qualHistogram[i];

		for (int i = 80; i <= 104; i++)
			phred64Total += qualHistogram[i];

		if (phred33Total == 0 && phred64Total > 0)
			return 64;

		if (phred64Total == 0 && phred33Total > 0)
			return 33;

		return 0;
	}

	public void open(File input) throws IOException {
		posTrackInputStream = new PositionTrackingInputStream(new FileInputStream(input), input.length());

		InputStream contentInputStream = CompressionFormat.wrapStreamForParsing(posTrackInputStream, input.getName());

		reader = new BufferedReader(new InputStreamReader(contentInputStream), 32768);

		if (phredOffset == 0) {
			deque.clear();
			qualHistogram = new int[256];
			long prereadsBytes = 0;

			for (int i = 0; i < PREREAD_COUNT; i++) {
				parseOne();
				if (current == null) {
					break; // EOF — stop preread early
				}
				deque.add(current);
				accumulateHistogram(current);
				prereadsBytes += current.getRecordLength();
				if (prereadsBytes >= PREREAD_MAX_BYTES) {
					break; // byte budget exhausted — enough data for phred detection
				}
			}
		}
		parseOne();
	}

	/**
	 * Open a FASTA file (optionally compressed) and prime the first record.
	 * Phred detection is skipped — callers must supply phredOffset=33 beforehand.
	 */
	public void openFasta(File input) throws IOException {
		posTrackInputStream = new PositionTrackingInputStream(new FileInputStream(input), input.length());
		InputStream contentInputStream = CompressionFormat.wrapStreamForParsing(posTrackInputStream, input.getName());
		reader = new BufferedReader(new InputStreamReader(contentInputStream), 32768);
		parseOne();
	}

	/**
	 * Return true if the first non-blank line of the (optionally compressed) file
	 * starts with '>' (FASTA format).
	 */
	public static boolean isFasta(File input) throws IOException {
		try (InputStream raw = new FileInputStream(input);
				InputStream decompressed = CompressionFormat.wrapStreamForParsing(raw, input.getName());
				BufferedReader br = new BufferedReader(new InputStreamReader(decompressed))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (!line.isEmpty())
					return line.charAt(0) == '>';
			}
		}
		return false;
	}

	public void close() throws IOException {
		// open() may not have completed (e.g. FileInputStream threw before reader
		// was assigned) — close() must still be safe to call in that case.
		if (reader != null)
			reader.close();
	}

	public boolean hasNext() {
		return (!deque.isEmpty()) || (current != null);
	}

	public FastqRecord next() throws IOException {
		if (deque.isEmpty()) {
			FastqRecord current = this.current;
			parseOne();

			return current;
		} else {
			FastqRecord rec = deque.poll();

			if (rec != null)
				rec.setPhredOffset(phredOffset);

			return rec;
		}
	}

}
