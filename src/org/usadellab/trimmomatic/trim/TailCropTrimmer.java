package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TailCropTrimmer extends AbstractSingleRecordTrimmer
{
	private int bases;
	private int maxLength=Integer.MAX_VALUE/2;

	public TailCropTrimmer(String args)
	{
		String arg[]=args.split(":");
	
		bases=Integer.parseInt(arg[0]);
	
		if(arg.length>1)
			maxLength=Integer.parseInt(arg[1]);
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
		int len=in.getSequence().length();
	
		int toTrim=bases;
		int overLen=len-toTrim-maxLength;
	
		if(overLen>0)
			toTrim+=overLen;

		if(len<=toTrim)
			return null;
	
		if(toTrim==0)
			return in;
	
		return new FastqRecord(in,0,len-toTrim);
	}

}
