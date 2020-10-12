package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class HeadCropTrimmer extends AbstractSingleRecordTrimmer
{
	private int pos;

	public HeadCropTrimmer(String args)
	{
		pos=Integer.parseInt(args);
	}

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int len=in.getSequence().length();
	
		if(len<=pos)
			return null;
	
		return new FastqRecord(in,pos,len-pos);
	}

}
