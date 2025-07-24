package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class ToPhred33Trimmer extends AbstractSingleRecordTrimmer
{

	public ToPhred33Trimmer(String args)
	{
	
	}

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		if(in.getPhredOffset()==33)
			return in;
	
		String sequence=in.getSequence();
		String quality=in.getQuality();
		
		StringBuilder newQuality=new StringBuilder();
		
		for(int i=0;i<quality.length();i++)
			{
			char newCh=(char)(quality.charAt(i)-31);
			newQuality.append(newCh);
			}
		
		return new FastqRecord(in, sequence, newQuality.toString(), 33);
	}

}
