package org.usadellab.trimmomatic.threading.trimlog;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public abstract class TrimLogCollector
{
	public static TrimLogCollector makeTrimLogCollector(boolean useWorker, int threads, File trimLog, ExceptionHolder exceptionHolder) throws IOException
	{
		if (trimLog==null)
			return null;
			
		PrintStream trimLogStream = new PrintStream(trimLog);
		
		if (useWorker)
			return new SelfThreadedTrimLogCollector(trimLogStream, threads * 5, exceptionHolder);
		else
			return new ParasiteTrimLogCollector(trimLogStream);
		
	}
	
	private AtomicBoolean complete;
	private PrintStream trimLogStream;
	
	TrimLogCollector(PrintStream trimLogStream)
	{
		this.trimLogStream = trimLogStream;	
		this.complete=new AtomicBoolean();
	}
	
	public boolean isComplete()
	{
		return complete.get();
	}
	
	protected void setCompleted()
	{
		complete.set(true);
	}

	protected void logRec(TrimLogRecord rec)
	{
		trimLogStream.printf("%s %d %d %d %d\n",rec.getReadName(),rec.getLength(),rec.getStartPos(),rec.getEndPos(),rec.getTrimTail());
	}
	
	public void close() throws Exception
	{
		trimLogStream.close();
	}
	
	public abstract void put(Future<BlockOfRecords> future) throws Exception;

	
}
