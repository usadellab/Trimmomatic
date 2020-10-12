package org.usadellab.trimmomatic.util;

import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.atomic.AtomicLong;

public class PositionTrackingInputStream extends InputStream
{
	private InputStream wrappedStream;
	private AtomicLong currentPosition;
	private AtomicLong markedPosition;

	public PositionTrackingInputStream(InputStream wrappedStream)
	{
		this.wrappedStream=wrappedStream;
		
		this.currentPosition=new AtomicLong();
		this.markedPosition=new AtomicLong();
	}
	
	@Override
	public int available() throws IOException
	{
		return wrappedStream.available();
	}

	@Override
	public void close() throws IOException
	{
		wrappedStream.close();
	}

	@Override
	public void mark(int readlimit)
	{
		wrappedStream.mark(readlimit);
		markedPosition.set(currentPosition.get());
	}

	@Override
	public boolean markSupported()
	{
		return wrappedStream.markSupported();
	}

	@Override
	public int read(byte[] arg0, int arg1, int arg2) throws IOException
	{
		int read=wrappedStream.read(arg0, arg1, arg2);
		
		if(read>0)
			currentPosition.addAndGet(read);
		
		return read;
	}

	@Override
	public int read(byte[] b) throws IOException
	{
		int read=wrappedStream.read(b);
		
		if(read>0)
			currentPosition.addAndGet(read);
	
		return read;	
	}

	@Override
	public void reset() throws IOException
	{
		wrappedStream.reset();
		
		currentPosition.set(markedPosition.get());		
	}

	@Override
	public long skip(long arg0) throws IOException
	{
		long read=wrappedStream.skip(arg0);
	
		if(read>0)
			currentPosition.addAndGet(read);

		return read;	
	}

	@Override
	public int read() throws IOException
	{
		int read=wrappedStream.read();
		
		if(read>=0)
			currentPosition.incrementAndGet();
		
		return read;
	}

	public long getPosition() 
	{
		return currentPosition.get();
	}
	
}
