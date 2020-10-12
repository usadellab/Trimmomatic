package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ParserWorker;
import org.usadellab.trimmomatic.threading.SerializerWorker;
import org.usadellab.trimmomatic.threading.TrimLogWorker;
import org.usadellab.trimmomatic.threading.TrimStatsWorker;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.trim.TrimmerFactory;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmomaticSE extends Trimmomatic
{

	/**
	 * Trimmomatic: The FASTQ trimmer
	 * 
	 * CROP:<LENGTH> Crop the read to specified length, by cutting off the right
	 * end LEADING:<QUAL> Trim the read by cutting off the left end while below
	 * specified quality TRAILING:<QUAL> Trim the read by cutting off the right
	 * end while below specified quality SLIDINGWINDOW:<QUAL>:<COUNT> Trim the
	 * read once the total quality of COUNT bases drops below QUAL, then trim
	 * trailing bases below QUAL
	 * 
	 * MINLEN:<LENGTH> Drop the read if less than specified length
	 */

	private Logger logger;

	public TrimmomaticSE(Logger logger)
	{
		this.logger=logger;
	}

	public void processSingleThreaded(FastqParser parser, FastqSerializer serializer, Trimmer trimmers[],
			PrintStream trimLogStream, PrintStream statsSummaryStream) throws IOException
	{
		TrimStats stats = new TrimStats();

		FastqRecord recs[] = new FastqRecord[1];
		FastqRecord originalRecs[] = new FastqRecord[1];

		while (parser.hasNext())
			{
			originalRecs[0] = recs[0] = parser.next();

			for (int i = 0; i < trimmers.length; i++)
				{
				try
					{
					recs = trimmers[i].processRecords(recs);
					}
				catch (RuntimeException e)
					{
					logger.errorln("Exception processing read: " + originalRecs[0].getName());
					throw e;
					}
				}

			if (recs[0] != null)
				{
				serializer.writeRecord(recs[0]);
				}

			stats.logPair(originalRecs, recs);

			if (trimLogStream != null)
				{
				for (int i = 0; i < originalRecs.length; i++)
					{
					int length = 0;
					int startPos = 0;
					int endPos = 0;
					int trimTail = 0;

					if (recs[i] != null)
						{
						length = recs[i].getSequence().length();
						startPos = recs[i].getHeadPos();
						endPos = length + startPos;
						trimTail = originalRecs[i].getSequence().length() - endPos;
						}

					trimLogStream.printf("%s %d %d %d %d\n", originalRecs[i].getName(), length, startPos, endPos,
							trimTail);
					}
				}
			}

		logger.infoln(stats.processStatsSE(statsSummaryStream));
	}

	public void processMultiThreaded(FastqParser parser, FastqSerializer serializer, Trimmer trimmers[],
			PrintStream trimLogStream, PrintStream statsSummaryStream, int threads) throws IOException
	{
		ArrayBlockingQueue<List<FastqRecord>> parserQueue = new ArrayBlockingQueue<List<FastqRecord>>(threads);
		ArrayBlockingQueue<Runnable> taskQueue = new ArrayBlockingQueue<Runnable>(threads * 2);
		ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads * 5);

		ParserWorker parserWorker = new ParserWorker(parser, parserQueue);
		Thread parserThread = new Thread(parserWorker);
		ThreadPoolExecutor taskExec = new ThreadPoolExecutor(threads, threads, 0, TimeUnit.SECONDS, taskQueue);
		SerializerWorker serializerWorker = new SerializerWorker(serializer, serializerQueue, 0);
		Thread serializerThread = new Thread(serializerWorker);

		ArrayBlockingQueue<Future<BlockOfRecords>> trimStatsQueue = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads * 5);
		TrimStatsWorker statsWorker = new TrimStatsWorker(trimStatsQueue);
		Thread statsThread = new Thread(statsWorker);

		ArrayBlockingQueue<Future<BlockOfRecords>> trimLogQueue = null;
		TrimLogWorker trimLogWorker = null;
		Thread trimLogThread = null;

		if (trimLogStream != null)
			{
			trimLogQueue = new ArrayBlockingQueue<Future<BlockOfRecords>>(threads * 5);
			trimLogWorker = new TrimLogWorker(trimLogStream, trimLogQueue);
			trimLogThread = new Thread(trimLogWorker);
			trimLogThread.start();
			}

		parserThread.start();
		serializerThread.start();
		statsThread.start();

		boolean done = false;

		List<FastqRecord> recs1 = null;

		try
			{
			while (!done)
				{
				recs1 = null;
				while (recs1 == null)
					recs1 = parserQueue.poll(1, TimeUnit.SECONDS);

				if (recs1 == null || recs1.size() == 0)
					done = true;

				BlockOfRecords bor = new BlockOfRecords(recs1, null);
				BlockOfWork work = new BlockOfWork(logger, trimmers, bor, false, trimLogStream != null);

				while (taskQueue.remainingCapacity() < 1)
					Thread.sleep(100);

				Future<BlockOfRecords> future = taskExec.submit(work);

				serializerQueue.put(future);
				trimStatsQueue.put(future);

				if (trimLogQueue != null)
					trimLogQueue.put(future);
				}

			parserThread.join();
			parser.close();

			taskExec.shutdown();
			taskExec.awaitTermination(1, TimeUnit.HOURS);

			serializerThread.join();
			if (trimLogThread != null)
				trimLogThread.join();

			statsThread.join();
			logger.infoln(statsWorker.getStats().processStatsSE(statsSummaryStream));
			}
		catch (InterruptedException e)
			{
			throw new RuntimeException(e);
			}

	}

	public void process(File input, File output, Trimmer trimmers[], int phredOffset, File trimLog, File statsSummaryFile, int threads)
			throws IOException
	{
		FastqParser parser = new FastqParser(phredOffset);
		parser.parse(input);

		if(phredOffset==0)
			{
			int phred=parser.determinePhredOffset();
			if(phred!=0)
				{
				logger.infoln("Quality encoding detected as phred"+phred);
				parser.setPhredOffset(phred);
				}
			else
				{
				logger.errorln("Error: Unable to detect quality encoding");
				System.exit(1);
				}
			}
		
		FastqSerializer serializer = new FastqSerializer();
		serializer.open(output);

		PrintStream trimLogStream = null;
		if (trimLog != null)
			trimLogStream = new PrintStream(trimLog);

		PrintStream statsSummaryStream = null;
		if(statsSummaryFile!=null)
			 statsSummaryStream = new PrintStream(statsSummaryFile);
		
		if (threads == 1)
			processSingleThreaded(parser, serializer, trimmers, trimLogStream, statsSummaryStream);
		else
			processMultiThreaded(parser, serializer, trimmers, trimLogStream, statsSummaryStream, threads);

		serializer.close();

		if (trimLogStream != null)
			trimLogStream.close();
		
		if(statsSummaryStream != null)
			statsSummaryStream.close();
	}

	public static boolean run(String[] args) throws IOException
	{
		int argIndex = 0;
		int phredOffset = 0;
		int threads = 0;

		boolean badOption = false;

		File trimLog = null;
		File statsSummary = null;
		boolean quiet=false;
		boolean showVersion=false;		

		List<String> nonOptionArgs=new ArrayList<String>();
		
		while (argIndex < args.length)
			{
			String arg = args[argIndex++];
			
			if(arg.startsWith("-"))
				{			
				if (arg.equals("-phred33"))
					phredOffset = 33;
				else if (arg.equals("-phred64"))
					phredOffset = 64;
				else if (arg.equals("-threads"))
					threads = Integer.parseInt(args[argIndex++]);
				else if (arg.equals("-trimlog"))
					{
					if (argIndex < args.length)
						trimLog = new File(args[argIndex++]);
					else
						badOption = true;
					}
				else if (arg.equals("-summary"))
					{
					if (argIndex < args.length)
						statsSummary = new File(args[argIndex++]);
					else
						badOption = true;
					}
				else if (arg.equals("-quiet"))
					quiet=true;
				else if (arg.equals("-version"))
					showVersion=true; 
				else
					{
					System.err.println("Unknown option " + arg);
					badOption = true;
					}
				}
			else
				nonOptionArgs.add(arg);
			}

		if(showVersion)
			Trimmomatic.showVersion();
		
		if ((nonOptionArgs.size() < 3) || badOption)
			return showVersion;

		Logger logger=new Logger(true,true,!quiet);
		
		logger.infoln("TrimmomaticSE: Started with arguments:");
		for (String arg : args)
			logger.info(" " + arg);
		logger.infoln();
		
		if(threads==0)
			{
			threads=calcAutoThreadCount();
			logger.infoln("Automatically using "+threads+" threads");
			}
		
		Iterator<String> nonOptionArgsIter=nonOptionArgs.iterator();
		
		File input = new File(nonOptionArgsIter.next());
		File output = new File(nonOptionArgsIter.next());

		Trimmer trimmers[]=createTrimmers(logger, nonOptionArgsIter);

		TrimmomaticSE tm = new TrimmomaticSE(logger);
		tm.process(input, output, trimmers, phredOffset, trimLog, statsSummary, threads);

		logger.infoln("TrimmomaticSE: Completed successfully");
		return true;
	}

	public static void main(String[] args) throws IOException
	{
		if(!run(args))
			{
			System.err
					.println("Usage: TrimmomaticSE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...");
			System.exit(1);
			}
	}

}
