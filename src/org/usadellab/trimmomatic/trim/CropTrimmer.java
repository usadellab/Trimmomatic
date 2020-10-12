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
        

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		if(in.getSequence().length()<len)
			return in;
	
		return new FastqRecord(in,0,len);
	}

}
