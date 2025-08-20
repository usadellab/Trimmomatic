package org.usadellab.trimmomatic.threading.parser;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class SelfThreadedParser extends Parser implements Runnable {
	private ArrayBlockingQueue<List<FastqRecord>> parserQueue;

	private ExceptionHolder exceptionHolder;
	private Thread thread;

	public SelfThreadedParser(FastqParser parser, int buffers, ExceptionHolder exceptionHolder) {
		super(parser);

		this.exceptionHolder = exceptionHolder;
		this.parserQueue = new ArrayBlockingQueue<List<FastqRecord>>(buffers);

		thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}

	public List<FastqRecord> poll() throws Exception {
		exceptionHolder.rethrow();

		return parserQueue.poll(100, TimeUnit.MILLISECONDS);
	}

	public void close() throws Exception {
		while (thread.isAlive()) {
			exceptionHolder.rethrow();
			thread.join(100);
		}

		super.close();
	}

	@Override
	public void run() {
		try {
			List<FastqRecord> recs = parseBlock();

			while (recs.size() > 0) {
				parserQueue.put(recs);
				recs = parseBlock();
			}

		} catch (Exception e) {
			Exception pe = new Exception("Parser Exception", e);
			exceptionHolder.setException(pe);
		} finally {
			setCompleted();

			try {
				parserQueue.put(new ArrayList<FastqRecord>());
			} catch (InterruptedException e) {
				exceptionHolder.setException(e);
			}
		}

	}

}
