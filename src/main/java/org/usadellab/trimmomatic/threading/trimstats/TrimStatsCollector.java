package org.usadellab.trimmomatic.threading.trimstats;

import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public abstract class TrimStatsCollector {
	public static TrimStatsCollector makeTrimStatsCollector(boolean useWorker, int threads,
			ExceptionHolder exceptionHolder) {
		if (useWorker)
			return new SelfThreadedTrimStatsCollector(threads, exceptionHolder);
		else
			return new ParasiteTrimStatsCollector();
	}

	private AtomicBoolean complete;

	TrimStatsCollector() {
		this.complete = new AtomicBoolean();
	}

	public boolean isComplete() {
		return complete.get();
	}

	protected void setCompleted() {
		complete.set(true);
	}

	public abstract void put(Future<BlockOfRecords> future) throws Exception;

	public abstract void close() throws Exception;

	public abstract TrimStats getStats();

}
