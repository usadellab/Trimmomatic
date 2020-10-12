package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MaxLenTrimmer extends AbstractSingleRecordTrimmer
{
	private int maxLen;

	public MaxLenTrimmer(String args)
	{
		maxLen=Integer.parseInt(args);
	}

        public MaxLenTrimmer(int maxLen) {
            this.maxLen = maxLen;
        }

        

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		if(in.getSequence().length()<=maxLen)
			return in;
		
		return null;
	}

}
