package org.usadellab.trimmomatic.threading;

import java.io.PrintStream;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

public class TrimLogWorker implements Runnable
{
	private PrintStream trimLogStream;
	private ArrayBlockingQueue<Future<BlockOfRecords>> logQueue;

	private AtomicBoolean complete;

	public TrimLogWorker(PrintStream trimLogStream, ArrayBlockingQueue<Future<BlockOfRecords>> logQueue)
	{
		this.trimLogStream = trimLogStream;
		this.logQueue=logQueue;
		this.complete=new AtomicBoolean();
	}
	
	public boolean isComplete()
	{
		return complete.get();
	}
	
	private void logRec(TrimLogRecord rec)
	{
		trimLogStream.printf("%s %d %d %d %d\n",rec.getReadName(),rec.getLength(),rec.getStartPos(),rec.getEndPos(),rec.getTrimTail());
	
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
			e.printStackTrace();
			throw new RuntimeException(e);
			}
		finally
			{
			complete.set(true);
			}
	}

}
