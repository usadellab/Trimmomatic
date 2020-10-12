package org.usadellab.trimmomatic.threading;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.TrimStats;

public class TrimStatsWorker implements Runnable
{
	private ArrayBlockingQueue<Future<BlockOfRecords>> logQueue;
	private AtomicBoolean complete;
	
	private TrimStats stats;

	public TrimStatsWorker(ArrayBlockingQueue<Future<BlockOfRecords>> logQueue)
	{
		this.logQueue=logQueue;
		this.complete=new AtomicBoolean();
		
		stats=new TrimStats();
	}
	
	public boolean isComplete()
	{
		return complete.get();
	}
		
	public TrimStats getStats()
	{
		return stats;
	}
	
	@Override
	public void run()
	{
		try
			{
			Future<BlockOfRecords> future=logQueue.take();
			BlockOfRecords bor=future.get();
			TrimStats st=bor.getStats();
			
			while(st!=null)
				{
				stats.merge(st);
				
				future=logQueue.take();
				bor=future.get();
				st=bor.getStats();
				}
				
			}
		catch (Exception e)
			{
			e.printStackTrace();
			throw new RuntimeException(e);
			}
		finally
			{
			complete.set(true);
			}
	}

}
