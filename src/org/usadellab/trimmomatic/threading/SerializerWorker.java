package org.usadellab.trimmomatic.threading;

import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;

public class SerializerWorker implements Runnable
{
	private FastqSerializer serializer;
	private ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue;
	private int recIndex;
	private AtomicBoolean complete;

	public SerializerWorker(FastqSerializer serializer, ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue, int recIndex)
	{
		this.serializer = serializer;
		this.serializerQueue = serializerQueue;
		this.recIndex = recIndex;
		
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
			Future<BlockOfRecords> future=serializerQueue.take();
			BlockOfRecords bor=future.get();
			List<FastqRecord> recs = bor.getTrimmedRecs().get(recIndex);

			while(recs!=null)
				{
				for(FastqRecord rec: recs)
					serializer.writeRecord(rec);
				
				future=serializerQueue.take();
				bor=future.get();
				recs = bor.getTrimmedRecs().get(recIndex);
				}
				
			}
		catch (Exception e)
			{
			e.printStackTrace();
			throw new RuntimeException(e);
			}
		finally
			{
			complete.set(true);
			}
	}

}
