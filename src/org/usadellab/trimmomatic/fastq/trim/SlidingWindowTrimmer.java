package org.usadellab.trimmomatic.fastq.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class SlidingWindowTrimmer extends AbstractSingleRecordTrimmer
{
	private int windowLength;
	private float requiredQuality;
	private float totalRequiredQuality;

	public SlidingWindowTrimmer(String args)
	{
		String arg[]=args.split(":");		
		windowLength=Integer.parseInt(arg[0]);
		requiredQuality=Float.parseFloat(arg[1]);		
		totalRequiredQuality=requiredQuality*windowLength; // Convert to total
	}

    public SlidingWindowTrimmer(int windowLength, float requiredQuality) {
        this.windowLength = windowLength;
        this.requiredQuality = requiredQuality;
        totalRequiredQuality=requiredQuality*windowLength; // Convert to total
    }

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int quals[]=in.getQualityAsInteger(true);
		
		if(quals.length<windowLength)
			return null;
		
		int total=0;
		for(int i=0;i<windowLength;i++)
			total+=quals[i];
		
		if(total<totalRequiredQuality)
			return null;
		
		int lengthToKeep=quals.length;
		
		for(int i=0;i<quals.length;i++)
			{
			int qual=0;
			if(i+windowLength<quals.length)
				qual=quals[i+windowLength];
			
			total=total-quals[i]+qual;
			if(total<totalRequiredQuality)
				{
				lengthToKeep=i+1;
				break;
				}
			}
		
		/*
		if(lengthToKeep==quals.length)
			{
			// Shrink window at end 
		
			float tmpTotalRequiredQuality=totalRequiredQuality;
		
			for(int i=quals.length-windowLength;i<quals.length;i++)
				{
				//total=total-quals[i];
				//tmpTotalRequiredQuality-=requiredQuality;
				
				total=total-quals[i]+quals[quals.length-1];
				tmpWindowLength--;
			
				if(total<tmpTotalRequiredQuality)
					{
					lengthToKeep=i+tmpWindowLength;
					break;
					}
				}
			}
		*/
		/*
		int i=lengthToKeep;
		
		int lastBaseQuality=quals[i-1];
		while(lastBaseQuality < requiredQuality && i > 1)
			{
			i--;
			lastBaseQuality=quals[i-1];
			}
		*/
		
		int i=lengthToKeep;
		/*
		while(i < quals.length)
			{
			int baseQuality=quals[i];
			if(baseQuality<requiredQuality)
				break;		
			i++;
			}
		*/
		
		if(i<1)
			return null;
		
		if(i<quals.length)
			return new FastqRecord(in,0,i);
		
		return in;
	}

}
