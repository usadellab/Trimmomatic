package org.usadellab.trimmomatic.threading.pipeline;

import java.util.concurrent.Future;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public abstract class Pipeline {
	public static Pipeline makePipeline(int threads, ExceptionHolder exceptionHolder) {
		if (threads == 1)
			return new ParasiteSerialPipeline();
		else
			return new ThreadedPipeline(threads, exceptionHolder);
	}

	public abstract Future<BlockOfRecords> submit(BlockOfWork work) throws Exception;

	public abstract void close() throws InterruptedException;

}
