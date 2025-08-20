package org.usadellab.trimmomatic.threading;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.trimlog.TrimLogRecord;

public class BlockOfRecords implements Future<BlockOfRecords> {
	private List<FastqRecord> originalRecs1;
	private List<FastqRecord> originalRecs2;

	private List<TrimLogRecord> trimLogRec;

	private TrimStats stats;

	public BlockOfRecords(List<FastqRecord> originalRecs1, List<FastqRecord> originalRecs2) {
		this.originalRecs1 = originalRecs1;
		this.originalRecs2 = originalRecs2;
	}

	public List<FastqRecord> takeOriginalRecs1() {
		List<FastqRecord> or1 = originalRecs1;
		originalRecs1 = null;
		return or1;
	}

	public List<FastqRecord> takeOriginalRecs2() {
		List<FastqRecord> or2 = originalRecs2;
		originalRecs2 = null;
		return or2;
	}

	public List<TrimLogRecord> getTrimLogRecs() {
		return trimLogRec;
	}

	public void setTrimLogRecs(List<TrimLogRecord> trimLogRec) {
		this.trimLogRec = trimLogRec;
	}

	public TrimStats getStats() {
		return stats;
	}

	public void setStats(TrimStats stats) {
		this.stats = stats;
	}

	@Override
	public boolean cancel(boolean mayInterruptIfRunning) {
		return false;
	}

	@Override
	public boolean isCancelled() {
		return false;
	}

	@Override
	public boolean isDone() {
		return true;
	}

	@Override
	public BlockOfRecords get() throws InterruptedException, ExecutionException {
		return this;
	}

	@Override
	public BlockOfRecords get(long timeout, TimeUnit unit)
			throws InterruptedException, ExecutionException, TimeoutException {
		return this;
	}

}
