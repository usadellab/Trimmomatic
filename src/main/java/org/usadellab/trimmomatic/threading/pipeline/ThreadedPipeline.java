package org.usadellab.trimmomatic.threading.pipeline;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class ThreadedPipeline extends Pipeline
{
	private ExceptionHolder exceptionHolder;

	ArrayBlockingQueue<Runnable> taskQueue;
	ThreadPoolExecutor taskExec;
	

	public ThreadedPipeline(int threads, ExceptionHolder exceptionHolder)
	{	
		this.exceptionHolder=exceptionHolder;
		
		taskQueue = new ArrayBlockingQueue<Runnable>(threads * 2);
		
		ThreadFactory tf = new ThreadFactory() {
			public Thread newThread(Runnable r) 
				{
				Thread t = Executors.defaultThreadFactory().newThread(r);
				t.setDaemon(true);
				return t;
				} 
			};

		taskExec = new ThreadPoolExecutor(threads, threads, 0, TimeUnit.SECONDS, taskQueue, tf);		
	}
	
	public Future<BlockOfRecords> submit(BlockOfWork work) throws Exception
	{
		while (taskQueue.remainingCapacity() < 1)
			{			
			Thread.sleep(100);
			exceptionHolder.rethrow();
			}

		Future<BlockOfRecords> future = taskExec.submit(work);
		return future;
	}
	
	public void close() throws InterruptedException
	{
		taskExec.shutdown();
		taskExec.awaitTermination(1, TimeUnit.HOURS);
	}

	
}
