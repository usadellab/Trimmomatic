package org.usadellab.trimmomatic.threading.parser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public abstract class Parser {
	public static final int BLOCK_MAX_RECORDS = 32768; // 32k records per block

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
		List<FastqRecord> recs = new ArrayList<FastqRecord>(BLOCK_MAX_RECORDS);

		while (parser.hasNext()) {
			FastqRecord rec = parser.next();
			recs.add(rec);

			if (recs.size() >= BLOCK_MAX_RECORDS)
				return recs;
		}

		return recs;
	}

	public abstract List<FastqRecord> poll() throws Exception; // Empty means EOF, null means not ready (yet)

}
