package org.usadellab.trimmomatic.threading;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class ParserWorker implements Runnable
{
	public static final int BLOCKSIZE=1000;

	private FastqParser parser;
	private ArrayBlockingQueue<List<FastqRecord>> parserQueue;
	private AtomicBoolean complete;
	
	public ParserWorker(FastqParser parser, ArrayBlockingQueue<List<FastqRecord>> parserQueue)
	{
		this.parser=parser;
		this.parserQueue=parserQueue;
		
		this.complete=new AtomicBoolean();
	}
	
	public boolean isComplete()
	{
		return complete.get();
	}
	
	@Override
	public void run() 
	{
		try
			{
			List<FastqRecord> recs=new ArrayList<FastqRecord>(BLOCKSIZE);
	
			while(parser.hasNext())
				{
				recs.add(parser.next());
				if(recs.size()>=BLOCKSIZE)
					{
					parserQueue.put(recs);
					recs=new ArrayList<FastqRecord>();
					}
				}
			
			if(recs.size()>0)
				parserQueue.put(recs);
			}
		catch(IOException e)
			{
			e.printStackTrace();
			throw new RuntimeException(e);
			}
		catch(InterruptedException e)
			{
			e.printStackTrace();
			throw new RuntimeException(e);
			}
		finally
			{
			complete.set(true);
			try
				{
				parserQueue.put(new ArrayList<FastqRecord>());
				}
			catch(InterruptedException e)
				{
				e.printStackTrace();
				throw new RuntimeException(e);
				} 
			}
	}
	
	
}
