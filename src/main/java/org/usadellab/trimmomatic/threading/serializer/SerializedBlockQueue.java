package org.usadellab.trimmomatic.threading.serializer;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.util.compression.ParallelCompressor;
import org.usadellab.trimmomatic.threading.ExceptionHolder;
import org.usadellab.trimmomatic.util.compression.BlockData;
import org.usadellab.trimmomatic.util.compression.GzipBlockData;
import org.usadellab.trimmomatic.util.compression.UncompressedBlockData;

public class SerializedBlockQueue
{	
	private ArrayBlockingQueue<SerializedBlock> compressionQueue;
	private ParallelCompressor compressor;
	
	private UncompressedBlockData previous;
	
	private ArrayBlockingQueue<SerializedBlock> outputQueue;

	public SerializedBlockQueue(int buffers, ParallelCompressor compressor)
	{
		this.compressor = compressor; // Can be null
	
		this.compressionQueue = new ArrayBlockingQueue<SerializedBlock>(buffers);		
		this.outputQueue = new ArrayBlockingQueue<SerializedBlock>(buffers);		
	}
	
	public ParallelCompressor getParallelCompressor()
	{
		return compressor;
	}
	
	// Called from main thread only
	public void put(SerializedBlock block, ExceptionHolder exceptionHolder) throws Exception
	{
		while(!compressionQueue.offer(block, 100, TimeUnit.MILLISECONDS))
			exceptionHolder.rethrow();
	
		while(!outputQueue.offer(block, 100, TimeUnit.MILLISECONDS))
			exceptionHolder.rethrow();	
	}
	
	// Called from pipeline threads (main or worker)
	public void pollCompressible() throws Exception
	{
		while (pollCompressibleOnce());
	}
	
	
	private boolean pollCompressibleOnce() throws Exception
	{
	
		UncompressedBlockData currentForCompression=null;
		UncompressedBlockData previousForCompression=null;

		SerializedBlock blockToCompress=null;
			
		boolean progress=false;
			
		synchronized(this)
			{
			SerializedBlock next = compressionQueue.peek();
												// 	Can do these checks initially before Synchronized
			if (next!=null && next.isCompressible())	
				{
				next = compressionQueue.poll();
				progress = true;
					
				if (compressor==null)
					next.setData(next.getUncompressedData());
				else
					{
					blockToCompress=next;						
					currentForCompression = next.getUncompressedData();	
					
					previousForCompression = previous;					
					previous = currentForCompression;
					
					compressor.updateChecksumPreCompression(currentForCompression);
					}
				}
			}
			
		if (blockToCompress!=null)
			{
			// 	Drop lock for actual compression part
			BlockData compressedData = compressor.compress(previousForCompression, currentForCompression);				
			blockToCompress.setData(compressedData);
			}
			
		return progress;
			
	}
	

	public SerializedBlock pollOutput()
	{
		SerializedBlock block = outputQueue.poll();
	  	
		if (block!=null && compressor!=null)
			{
			try
				{
				BlockData blockData = block.get();
				if (blockData!=null)
					compressor.updateChecksumPostCompression(blockData);
				}
			catch (Exception e) 
				{ 
				throw new RuntimeException(e); 
				}
			}
		
	  	return block;
	}
	
	public SerializedBlock takeOutput() throws InterruptedException
	{
		SerializedBlock block = outputQueue.take();

		if (compressor!=null)
			{
			try
				{
				BlockData blockData = block.get();
				if (blockData!=null)
					compressor.updateChecksumPostCompression(blockData);
				}
			catch (Exception e) 
				{ 
				throw new RuntimeException(e); 
				}
			}
		
		return block;
	}
	
	
}
