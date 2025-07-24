package org.usadellab.trimmomatic.trim;

import java.util.BitSet;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BaseCountTrimmer extends AbstractSingleRecordTrimmer
{
	private int minCount=0;
	private Integer maxCount=null;
	private BitSet baseSet;

	public BaseCountTrimmer(String args)
	{
		String split[]=args.split(":");
		baseSet=new BitSet();
		
		String bases=split[0];
		
		for(int i=0;i<bases.length();i++)
			{
			char c=bases.charAt(i);
			baseSet.set(c);
			}	
		
		if(split.length>1)
			{		
			minCount=Integer.parseInt(split[1]);
			
			if(split.length>2)
				maxCount=new Integer(split[2]);
			}
	}
		
       
	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int count=0;
		
		String seq=in.getSequence();
		
		for(int i=0;i<seq.length();i++)
			{
			char c=seq.charAt(i);
			if(baseSet.get(c))
				count++;
			}			
	
		if(count<minCount)
			return null;
		
		if(maxCount!=null && count>maxCount)
			return null;
		
		return in;
	}

}
