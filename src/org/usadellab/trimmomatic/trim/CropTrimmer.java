package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class CropTrimmer extends AbstractSingleRecordTrimmer
{
	private int len;

	public CropTrimmer(String args)
	{
		len=Integer.parseInt(args);
	}

        public CropTrimmer(int len) {
            this.len = len;
        }        

        
        
/*
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
 */
        
        
	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		if(in.getSequence().length()<len)
			return in;
	
		return new FastqRecord(in,0,len);
	}

}
