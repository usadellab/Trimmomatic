package org.usadellab.trimmomatic.threading.trimstats;

import java.io.IOException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.threading.BlockOfRecords;

public class ParasiteTrimStatsCollector extends TrimStatsCollector
{
	private TrimStats stats;

	public ParasiteTrimStatsCollector()
	{
		stats=new TrimStats();
	}

	public void put(Future<BlockOfRecords> future) throws Exception
	{
		BlockOfRecords bor = future.get();
		TrimStats recs = bor.getStats();
		if (recs!=null)
			stats.merge(recs);
	}
	
	public void close() throws InterruptedException
	{		
	}
	
	public TrimStats getStats()
	{
		return stats;
	}

}
