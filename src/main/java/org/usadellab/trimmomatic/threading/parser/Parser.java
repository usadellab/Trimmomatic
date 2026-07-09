package org.usadellab.trimmomatic.threading.parser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public abstract class Parser implements AutoCloseable {
	public static final int  BLOCK_MAX_RECORDS = 32768;             // 32k records per block
	public static final long BLOCK_MAX_BYTES   = 64L * 1024 * 1024; // 64 MB per block (long-read safety)

	public static Parser makeParser(boolean useWorker, int buffers, FastqParser rawParser,
			ExceptionHolder exceptionHolder) {
		if (useWorker)
			return new SelfThreadedParser(rawParser, buffers, exceptionHolder);
		else
			return new ParasiteSerialParser(rawParser);
	}

	private FastqParser parser;
	private AtomicBoolean complete;

	Parser(FastqParser parser) {
		this.parser = parser;
		this.complete = new AtomicBoolean();
	}

	public boolean isComplete() {
		return complete.get();
	}

	protected void setCompleted() {
		complete.set(true);
	}

	public void close() throws Exception {
		parser.close();
	}

	protected List<FastqRecord> parseBlock() throws IOException {
		// Initial capacity is capped at 1024 to avoid over-allocating for long-read
		// files where a single block may contain very few but very large records.
		List<FastqRecord> recs = new ArrayList<>(Math.min(BLOCK_MAX_RECORDS, 1024));
		long blockBytes = 0;

		while (parser.hasNext()) {
			FastqRecord rec = parser.next();
			recs.add(rec);
			blockBytes += rec.getRecordLength();

			if (recs.size() >= BLOCK_MAX_RECORDS || blockBytes >= BLOCK_MAX_BYTES) {
				return recs;
			}
		}

		return recs;
	}

	public abstract List<FastqRecord> poll() throws Exception; // Empty means EOF, null means not ready (yet)

}
