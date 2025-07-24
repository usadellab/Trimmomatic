package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TrailingTrimmer extends AbstractSingleRecordTrimmer
{
	private int qual;

	public TrailingTrimmer(String args)
	{
		qual=Integer.parseInt(args);
	}

        public TrailingTrimmer(int qual) {
            this.qual = qual;
        }        
        

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int quals[]=in.getQualityAsInteger(true);
		
		for(int i=quals.length-1;i>0;i--)
			{
			if(quals[i]>=qual)
				return new FastqRecord(in,0,i+1);
			}
		
		return null;
	}

}
