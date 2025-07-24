package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public abstract class AbstractSingleRecordTrimmer implements Trimmer
{

	@Override
	public FastqRecord[] processRecords(FastqRecord[] in)
	{
		if(in==null)
			return null; 

		FastqRecord out[]=new FastqRecord[in.length];
		
		for(int i=0;i<in.length;i++)
			{
			if(in[i]!=null)
				out[i]=processRecord(in[i]);
			}
	
		return out;
	}

	public abstract FastqRecord processRecord(FastqRecord in);
	
	
}
