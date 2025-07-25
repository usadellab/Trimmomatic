package org.usadellab.trimmomatic.threading.trimstats;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class SelfThreadedTrimStatsCollector extends TrimStatsCollector implements Runnable
{
	private ArrayBlockingQueue<Future<BlockOfRecords>> statsQueue;
	
	private TrimStats stats;

	private ExceptionHolder exceptionHolder;
	private Thread thread;
	
	public SelfThreadedTrimStatsCollector(int buffers, ExceptionHolder exceptionHolder)
	{
		this.exceptionHolder = exceptionHolder;
		this.statsQueue=new ArrayBlockingQueue<Future<BlockOfRecords>>(buffers);
		
		stats=new TrimStats();
		
		thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}
	
	public void put(Future<BlockOfRecords> future) throws Exception
	{
		while(!statsQueue.offer(future, 100, TimeUnit.MILLISECONDS))
			exceptionHolder.rethrow();
	}
	
	public void close() throws Exception
	{
		while (thread.isAlive())
			{
			exceptionHolder.rethrow();
			thread.join(100);
			}		
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
			Future<BlockOfRecords> future=statsQueue.take();
			BlockOfRecords bor=future.get();
			TrimStats st=bor.getStats();
			
			while(st!=null)
				{
				stats.merge(st);
				
				future=statsQueue.take();
				bor=future.get();
				st=bor.getStats();
				}
				
			}
		catch (Exception e)
			{
			Exception pe = new Exception("Trim Stats Collector Exception", e);			
			exceptionHolder.setException(pe);
			}
		finally
			{
			setCompleted();
			}
	}

}
