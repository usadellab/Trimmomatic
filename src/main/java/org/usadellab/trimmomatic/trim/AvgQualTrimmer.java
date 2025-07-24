package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class AvgQualTrimmer extends AbstractSingleRecordTrimmer
{
    private int qual;

    public AvgQualTrimmer(String args)
    {
            qual=Integer.parseInt(args);
    }

    public AvgQualTrimmer(int qual) {
        this.qual = qual;
    }
        
        

    @Override
    public FastqRecord processRecord(FastqRecord in)
    {
            String seq=in.getSequence();
            int quals[]=in.getQualityAsInteger(true);

            int total=0;
            
            for(int i=0;i<seq.length();i++)
                total+=quals[i];

            if(total<qual*seq.length())
              return null;
            
            return in;
    }

}
