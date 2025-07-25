package org.usadellab.trimmomatic.threading.serializer;

import java.io.IOException;
import java.io.OutputStream;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.util.compression.BlockData;

public class ParasiteSerializer extends Serializer
{
	private AtomicBoolean complete;

	public ParasiteSerializer(OutputStream stream, SerializedBlockQueue sbq) throws IOException
	{
		super(stream, sbq);
		
		this.complete=new AtomicBoolean();
	}
		
	public void pollWritable() throws Exception
	{
		SerializedBlock block = sbq.pollOutput();
		
		if (block==null)
			throw new RuntimeException("Missing expected SerializedBlock");
	
		if (!block.isDone())
			throw new RuntimeException("Incomplete SerializedBlock for ParasiteSerializer");
		
		BlockData blockData = block.get();						
		writeBuffer(blockData);			
	}

	public boolean isComplete()
	{
		return complete.get();
	}
		
	private void writeBuffer(BlockData blockData) throws IOException
	{
		if (blockData!=null)
			{
			stream.writeBlock(blockData);
			}
		else
			{
			complete.set(true);
			}
	}

}
