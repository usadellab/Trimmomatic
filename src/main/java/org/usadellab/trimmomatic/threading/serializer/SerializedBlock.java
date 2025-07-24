package org.usadellab.trimmomatic.threading.serializer;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.util.compression.BlockData;
import org.usadellab.trimmomatic.util.compression.UncompressedBlockData;

public class SerializedBlock implements Future<BlockData>
{
	boolean last;

	AtomicBoolean compressible;
	AtomicBoolean done;
	
	UncompressedBlockData ucData;
	BlockData data;
	
	public SerializedBlock(boolean last)
	{	
		this.last = last;
		
		compressible=new AtomicBoolean();
		done=new AtomicBoolean();
	}
		
	public boolean isLast()
	{
		return last;
	}
	
	public void setUncompressedData(UncompressedBlockData ucData)
	{
		synchronized(this)
			{
			this.ucData = ucData;
			compressible.set(true);
			notifyAll();
			}
	}
	
	public UncompressedBlockData getUncompressedData()
	{
		return ucData;
	}
		
	public boolean isCompressible()
	{
		return compressible.get();
	}
	
	
	public void setData(BlockData data)	
	{
		synchronized(this)
			{
			this.data=data;
			done.set(true);
			notifyAll();
			}
	}
		
	@Override
	public boolean cancel(boolean mayInterruptIfRunning)
	{
		return false;
	}

	@Override
	public boolean isCancelled()
	{
		return false;
	}

	@Override
	public boolean isDone()
	{
		return done.get();
	}

	@Override
	public BlockData get() throws InterruptedException, ExecutionException
	{
		synchronized (this)
			{
			while (!isDone())
				wait();
			
			return data;
			}
	}

	@Override
	public BlockData get(long timeout, TimeUnit unit) throws InterruptedException, ExecutionException, TimeoutException
	{
		return null;
	}
	
	
}
