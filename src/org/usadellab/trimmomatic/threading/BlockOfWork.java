package org.usadellab.trimmomatic.threading;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

public class BlockOfWork implements Callable<BlockOfRecords>
{
	private Logger logger;
	
	private Trimmer trimmers[];
	private BlockOfRecords bor;
	private boolean pe;
	private boolean trimLog;

	
	public BlockOfWork(Logger logger, Trimmer trimmers[], BlockOfRecords bor, boolean pe, boolean trimLog)
	{
		this.logger=logger;
		
		this.trimmers = trimmers;
		this.bor = bor;
		
		this.pe = pe;
		this.trimLog=trimLog;
	}

	
	private TrimLogRecord makeTrimLogRec(FastqRecord rec, FastqRecord originalRec)
	{
		int length=0;
		int startPos=0;
		int endPos=0;
		int trimTail=0;

		if(rec!=null)
			{
			length=rec.getSequence().length();
			startPos=rec.getHeadPos();
			endPos=length+startPos;
			trimTail=originalRec.getSequence().length()-endPos;
			}
		
		return new TrimLogRecord(originalRec.getName(), length, startPos, endPos, trimTail);
	}
	
	@Override
	public BlockOfRecords call() throws Exception
	{
		TrimStats stats=new TrimStats();
	
		if (pe)
			{
			List<FastqRecord> originalRecs1 = bor.getOriginalRecs1();
			List<FastqRecord> originalRecs2 = bor.getOriginalRecs2();

			int len1 = originalRecs1.size();
			int len2 = originalRecs2.size();

			if(len1 == 0 && len2==0)
				{
				List<List<FastqRecord>> trimmedRecs=new ArrayList<List<FastqRecord>>();
				
				for(int i=0;i<4;i++)
					trimmedRecs.add(null);
				
				bor.setTrimmedRecs(trimmedRecs, null);
				return bor;
				}
			
			int len = len1 < len2 ? len1 : len2;

			FastqRecord originalRecs[]=new FastqRecord[2];
			
			List<FastqRecord>trimmedRecs1P=new ArrayList<FastqRecord>();
			List<FastqRecord>trimmedRecs1U=new ArrayList<FastqRecord>();
			List<FastqRecord>trimmedRecs2P=new ArrayList<FastqRecord>();
			List<FastqRecord>trimmedRecs2U=new ArrayList<FastqRecord>();
			
			List<TrimLogRecord>trimLogList=null;
			if(trimLog)
				trimLogList=new ArrayList<TrimLogRecord>();
			
			for (int i = 0; i < len; i++)
				{
				originalRecs[0]=originalRecs1.get(i);
				originalRecs[1]=originalRecs2.get(i);
				FastqRecord recs[]=originalRecs;
				
				for(int j=0;j<trimmers.length;j++)
					{
					try
						{
						recs=trimmers[j].processRecords(recs);
						}
					catch (RuntimeException e)
						{
						logger.errorln("Exception processing reads: "+originalRecs[0].getName()+" and "+originalRecs[1].getName());
						throw e;
						}
					}
				
				if(recs[0]!=null && recs[1]!=null)
					{
					trimmedRecs1P.add(recs[0]);
					trimmedRecs2P.add(recs[1]);
					}
				else if(recs[0]!=null)
					trimmedRecs1U.add(recs[0]);
				else if(recs[1]!=null)
					trimmedRecs2U.add(recs[1]);
				
				stats.logPair(originalRecs, recs);
				
				if(trimLog)
					{
					if(originalRecs[0]!=null)
						trimLogList.add(makeTrimLogRec(recs[0], originalRecs[0]));

					if(originalRecs[1]!=null)
						trimLogList.add(makeTrimLogRec(recs[1], originalRecs[1]));
					}
				}
			
			List<List<FastqRecord>> trimmedRecsList=new ArrayList<List<FastqRecord>>();
			trimmedRecsList.add(trimmedRecs1P);
			trimmedRecsList.add(trimmedRecs1U);
			trimmedRecsList.add(trimmedRecs2P);
			trimmedRecsList.add(trimmedRecs2U);
			
			bor.setTrimmedRecs(trimmedRecsList, trimLogList);
			bor.setStats(stats);			
			}
		else
			{
			List<FastqRecord> originalRecsL = bor.getOriginalRecs1();

			int len = originalRecsL.size();

			if(len == 0)
				{
				List<List<FastqRecord>> trimmedRecs=new ArrayList<List<FastqRecord>>();
				for(int i=0;i<4;i++)
					trimmedRecs.add(null);
				
				bor.setTrimmedRecs(trimmedRecs, null);
				return bor;
				}
			
			FastqRecord originalRecs[]=new FastqRecord[1];
			
			List<FastqRecord>trimmedRecs=new ArrayList<FastqRecord>();
			
			List<TrimLogRecord>trimLogList=null;
			if(trimLog)
				trimLogList=new ArrayList<TrimLogRecord>();
			
			for (int i = 0; i < len; i++)
				{
				originalRecs[0]=originalRecsL.get(i);
				FastqRecord recs[]=originalRecs;
				
				for(int j=0;j<trimmers.length;j++)
					{
					try
						{
						recs=trimmers[j].processRecords(recs);
						}
					catch (RuntimeException e)
						{
						logger.errorln("Exception processing read: "+originalRecs[0].getName());
						e.printStackTrace();
						throw e;
						}
					}
				
				if(recs[0]!=null)
					trimmedRecs.add(recs[0]);

				stats.logPair(originalRecs, recs);
				
				if(trimLog)
					{
					if(originalRecs[0]!=null)
						trimLogList.add(makeTrimLogRec(recs[0], originalRecs[0]));
					}
				}
			
			List<List<FastqRecord>> trimmedRecsList=new ArrayList<List<FastqRecord>>();
			trimmedRecsList.add(trimmedRecs);
			
			bor.setTrimmedRecs(trimmedRecsList, trimLogList);
			bor.setStats(stats);
			}

		return bor;
	}
}
