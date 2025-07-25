package org.usadellab.trimmomatic.threading.pipeline;

import java.util.concurrent.Future;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;

public class ParasiteSerialPipeline extends Pipeline
{
	public ParasiteSerialPipeline()
	{
	
	}
	
	public Future<BlockOfRecords> submit(BlockOfWork work) throws Exception 
	{
		BlockOfRecords bor = work.process();
		return bor;
	}
	
	public void close()
	{		
	}
}
