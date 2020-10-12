package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TrimStats
{
	private static DecimalFormat formatter=new DecimalFormat("0.00");

	private long readsInput;
	private long readsSurvivingBoth;
	private long readsSurvivingForward;
	private long readsSurvivingReverse;
	
	public TrimStats()
	{
		readsInput=0;
		readsSurvivingBoth=0;
		readsSurvivingForward=0;
		readsSurvivingReverse=0;
	}
	
	public void merge(TrimStats other)
	{
		readsInput+=other.readsInput;
		readsSurvivingBoth+=other.readsSurvivingBoth;
		readsSurvivingForward+=other.readsSurvivingForward;
		readsSurvivingReverse+=other.readsSurvivingReverse;
	}
	
	
	public void logPair(FastqRecord originalRecs[], FastqRecord recs[])
	{
		if(originalRecs.length==1)
			{
			if(originalRecs[0]!=null)
				{
				readsInput++;
			
				if(recs[0]!=null)
					readsSurvivingForward++;
				}
			}
		else
			{
			if(originalRecs[0]!=null && originalRecs[1]!=null)
				{
				readsInput++;
			
				if(recs[0]!=null)
					{
					if(recs[1]!=null)
						readsSurvivingBoth++;
					else
						readsSurvivingForward++;
					}
				else if(recs[1]!=null)
					readsSurvivingReverse++;
				}
			
			}
	}
	
	public String processStatsSE(PrintStream statsSummaryStream)
	{
		long dropped=readsInput-readsSurvivingForward;
		
		double survivingForwardPercent=(100.0*readsSurvivingForward)/readsInput;
		double droppedPercent=(100.0*dropped)/readsInput;
		
		if(statsSummaryStream!=null)
			{
			statsSummaryStream.println("Input Reads: "+readsInput);
			statsSummaryStream.println("Surviving Reads: "+readsSurvivingForward);
			statsSummaryStream.println("Surviving Read Percent: "+formatter.format(survivingForwardPercent));
			statsSummaryStream.println("Dropped Reads: "+dropped);
			statsSummaryStream.println("Dropped Read Percent: "+formatter.format(droppedPercent));				
			}
				
		return "Input Reads: "+readsInput+
			   " Surviving: "+readsSurvivingForward+" ("+formatter.format(survivingForwardPercent)+
			   "%) Dropped: "+dropped+" ("+formatter.format(droppedPercent)+"%)";
	}
	
	public String processStatsPE(PrintStream statsSummaryStream)
	{
		long dropped=readsInput-readsSurvivingBoth-readsSurvivingForward-readsSurvivingReverse;
	
		double survivingBothPercent=(100.0*readsSurvivingBoth)/readsInput;
		double survivingForwardPercent=(100.0*readsSurvivingForward)/readsInput;
		double survivingReversePercent=(100.0*readsSurvivingReverse)/readsInput;
		double droppedPercent=(100.0*dropped)/readsInput;

		if(statsSummaryStream!=null)
			{
			statsSummaryStream.println("Input Read Pairs: "+readsInput);
			statsSummaryStream.println("Both Surviving Reads: "+readsSurvivingBoth);
			statsSummaryStream.println("Both Surviving Read Percent: "+formatter.format(survivingBothPercent));
			statsSummaryStream.println("Forward Only Surviving Reads: "+readsSurvivingForward);
			statsSummaryStream.println("Forward Only Surviving Read Percent: "+formatter.format(survivingForwardPercent));
			statsSummaryStream.println("Reverse Only Surviving Reads: "+readsSurvivingReverse);
			statsSummaryStream.println("Reverse Only Surviving Read Percent: "+formatter.format(survivingReversePercent));				
			statsSummaryStream.println("Dropped Reads: "+dropped);
			statsSummaryStream.println("Dropped Read Percent: "+formatter.format(droppedPercent));
			}
		
		return "Input Read Pairs: "+readsInput+
		   " Both Surviving: "+readsSurvivingBoth+" ("+formatter.format(survivingBothPercent)+
		   "%) Forward Only Surviving: "+readsSurvivingForward+" ("+formatter.format(survivingForwardPercent)+
		   "%) Reverse Only Surviving: "+readsSurvivingReverse+" ("+formatter.format(survivingReversePercent)+
		   "%) Dropped: "+dropped+" ("+formatter.format(droppedPercent)+"%)";

	}
}
