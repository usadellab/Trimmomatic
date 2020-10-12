package org.usadellab.trimmomatic;

import java.text.DecimalFormat;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TrimStats
{
	private static DecimalFormat formatter=new DecimalFormat("0.00");

	private int input;
	private int survivingBoth;
	private int survivingForward;
	private int survivingReverse;
	
	public TrimStats()
	{
		input=0;
		survivingBoth=0;
		survivingForward=0;
		survivingReverse=0;
	}
	
	public void merge(TrimStats other)
	{
		input+=other.input;
		survivingBoth+=other.survivingBoth;
		survivingForward+=other.survivingForward;
		survivingReverse+=other.survivingReverse;
	}
	
	
	public void logPair(FastqRecord originalRecs[], FastqRecord recs[])
	{
		if(originalRecs.length==1)
			{
			if(originalRecs[0]!=null)
				{
				input++;
			
				if(recs[0]!=null)
					survivingForward++;
				}
			}
		else
			{
			if(originalRecs[0]!=null && originalRecs[1]!=null)
				{
				input++;
			
				if(recs[0]!=null)
					{
					if(recs[1]!=null)
						survivingBoth++;
					else
						survivingForward++;
					}
				else if(recs[1]!=null)
					survivingReverse++;
				}
			
			}
	}
	
	public String getStatsSE()
	{
		int dropped=input-survivingForward;
		
		double survivingForwardPercent=(100.0*survivingForward)/input;
		double droppedPercent=(100.0*dropped)/input;
		
		return "Input Reads: "+input+
			   " Surviving: "+survivingForward+" ("+formatter.format(survivingForwardPercent)+
			   "%) Dropped: "+dropped+" ("+formatter.format(droppedPercent)+"%)";
	}
	
	public String getStatsPE()
	{
		int dropped=input-survivingBoth-survivingForward-survivingReverse;
	
		double survivingBothPercent=(100.0*survivingBoth)/input;
		double survivingForwardPercent=(100.0*survivingForward)/input;
		double survivingReversePercent=(100.0*survivingReverse)/input;
		double droppedPercent=(100.0*dropped)/input;

		return "Input Read Pairs: "+input+
		   " Both Surviving: "+survivingBoth+" ("+formatter.format(survivingBothPercent)+
		   "%) Forward Only Surviving: "+survivingForward+" ("+formatter.format(survivingForwardPercent)+
		   "%) Reverse Only Surviving: "+survivingReverse+" ("+formatter.format(survivingReversePercent)+
		   "%) Dropped: "+dropped+" ("+formatter.format(droppedPercent)+"%)";

	}
}
