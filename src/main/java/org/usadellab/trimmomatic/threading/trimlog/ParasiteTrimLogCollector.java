package org.usadellab.trimmomatic.threading.trimlog;

import java.io.PrintStream;
import java.util.List;
import java.util.concurrent.Future;

import org.usadellab.trimmomatic.threading.BlockOfRecords;

public class ParasiteTrimLogCollector extends TrimLogCollector {
	public ParasiteTrimLogCollector(PrintStream trimLogStream) {
		super(trimLogStream);
	}

	public void put(Future<BlockOfRecords> future) throws Exception {
		BlockOfRecords bor = future.get();
		List<TrimLogRecord> recs = bor.getTrimLogRecs();

		if (recs != null)
			logBlock(recs);
	}

	private void logBlock(List<TrimLogRecord> recs) {
		for (TrimLogRecord rec : recs)
			logRec(rec);
	}

}
