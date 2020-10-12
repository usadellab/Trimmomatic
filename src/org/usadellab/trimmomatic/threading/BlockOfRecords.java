package org.usadellab.trimmomatic.threading;

import java.util.List;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BlockOfRecords
{
	private List<FastqRecord>originalRecs1;
	private List<FastqRecord>originalRecs2;
	
	private List<List<FastqRecord>> trimmedRecs;

	private List<TrimLogRecord> trimLogRec;
	
	private TrimStats stats;
	
	public BlockOfRecords(List<FastqRecord> originalRecs1, List<FastqRecord> originalRecs2)
	{
		this.originalRecs1=originalRecs1;
		this.originalRecs2=originalRecs2;
	}

	public synchronized List<List<FastqRecord>> getTrimmedRecs()
	{
		return trimmedRecs;
	}

	public synchronized List<TrimLogRecord> getTrimLogRecs()
	{
		return trimLogRec;
	}
	
	public synchronized void setTrimmedRecs(List<List<FastqRecord>> trimmedRecs, List<TrimLogRecord> trimLogRec)
	{
		this.trimmedRecs = trimmedRecs;
		this.trimLogRec = trimLogRec;
	}

	public synchronized List<FastqRecord> getOriginalRecs1()
	{
		return originalRecs1;
	}

	public synchronized List<FastqRecord> getOriginalRecs2()
	{
		return originalRecs2;
	}

	public synchronized TrimStats getStats()
	{
		return stats;
	}

	public synchronized void setStats(TrimStats stats)
	{
		this.stats = stats;
	}
	
}
