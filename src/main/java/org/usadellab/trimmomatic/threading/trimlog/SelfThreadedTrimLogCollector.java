package org.usadellab.trimmomatic.threading.trimlog;

import java.io.PrintStream;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class SelfThreadedTrimLogCollector extends TrimLogCollector implements Runnable
{
	private ArrayBlockingQueue<Future<BlockOfRecords>> logQueue;

	private ExceptionHolder exceptionHolder;
	private Thread thread;
	
	public SelfThreadedTrimLogCollector(PrintStream trimLogStream, int buffers, ExceptionHolder exceptionHolder)
	{
		super(trimLogStream);
		this.logQueue=new ArrayBlockingQueue<Future<BlockOfRecords>>(buffers);
		
		this.exceptionHolder=exceptionHolder;
		
		thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}
	
	public void put(Future<BlockOfRecords> future) throws Exception
	{
		while(!logQueue.offer(future, 100, TimeUnit.MILLISECONDS))
			exceptionHolder.rethrow();
	}
	
	public void close() throws Exception
	{	
		while (thread.isAlive())
			{
			exceptionHolder.rethrow();
			thread.join(100);
			}
	
		super.close();
	}
	
	
	@Override
	public void run()
	{
		try
			{
			Future<BlockOfRecords> future=logQueue.take();
			BlockOfRecords bor=future.get();
			List<TrimLogRecord> recs = bor.getTrimLogRecs();

			while(recs!=null)
				{
				for(TrimLogRecord rec: recs)
					logRec(rec);
				
				future=logQueue.take();
				bor=future.get();
				recs = bor.getTrimLogRecs();
				}
				
			}
		catch (Exception e)
			{
			Exception pe = new Exception("Trim Log Collector Exception", e);			
			exceptionHolder.setException(pe);
			}
		finally
			{
			setCompleted();
			}
	}

}
