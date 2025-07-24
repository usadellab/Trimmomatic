package org.usadellab.trimmomatic.threading.serializer;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.threading.ExceptionHolder;
import org.usadellab.trimmomatic.util.Logger;
import org.usadellab.trimmomatic.util.compression.BlockOutputStream;
import org.usadellab.trimmomatic.util.compression.CompressionFormat;
import org.usadellab.trimmomatic.util.compression.ParallelCompressor;

public abstract class Serializer
{
	public static Serializer makeSerializer(Logger logger, boolean useWorker, boolean useParallelCompressor, Integer compressLevel, int buffers, File output, ExceptionHolder exceptionHolder) throws IOException
	{
		OutputStream stream = null;	
		ParallelCompressor bc = null;
			
		if (useParallelCompressor)
			{
			bc = CompressionFormat.parallelCompressorForSerializing(logger, output.getName(), compressLevel);
			stream = new FileOutputStream(output);
			}
		else
			{
			stream = CompressionFormat.wrapStreamForSerializing(new FileOutputStream(output), output.getName(), compressLevel);
			}
		
		SerializedBlockQueue sbq = new SerializedBlockQueue(buffers, bc);
		
		if(useWorker)
			return new SelfThreadedSerializer(stream, sbq, exceptionHolder);
		else
			return new ParasiteSerializer(stream, sbq);
	}

	
	protected BlockOutputStream stream;
	protected SerializedBlockQueue sbq;
	
	private AtomicBoolean complete;
	
	Serializer(OutputStream stream, SerializedBlockQueue sbq) throws IOException
	{
		this.sbq = sbq;
		
		ParallelCompressor bc = sbq.getParallelCompressor();
		if (bc!=null)
			this.stream=bc.wrapAndWriteHeader(stream);
		else
			this.stream=new BlockOutputStream(stream);
		
		this.sbq = sbq;
		
		this.complete=new AtomicBoolean();
	}
	
	public ParallelCompressor getParallelCompressor()
	{
		return sbq.getParallelCompressor();
	}
	
	public boolean isComplete()
	{
		return complete.get();
	}
	
	protected void setCompleted()
	{
		complete.set(true);
	}
	
	
	public void close() throws Exception
	{	
		ParallelCompressor bc = sbq.getParallelCompressor();
		if (bc!=null)
			bc.writeTrailer(stream);
	
		stream.close();
	}
	
	public void queueForWrite(SerializedBlock block, ExceptionHolder exceptionHolder) throws Exception
	{
		sbq.put(block, exceptionHolder);
	}
	
	public void pollCompressible() throws Exception
	{
		sbq.pollCompressible();
	}
	
	public abstract void pollWritable() throws Exception;
	
}
	
