package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Future;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.PairingValidator;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ExceptionHolder;
import org.usadellab.trimmomatic.threading.parser.Parser;
import org.usadellab.trimmomatic.threading.pipeline.Pipeline;
import org.usadellab.trimmomatic.threading.serializer.SerializedBlock;
import org.usadellab.trimmomatic.threading.serializer.Serializer;
import org.usadellab.trimmomatic.threading.trimlog.TrimLogCollector;
import org.usadellab.trimmomatic.threading.trimstats.TrimStatsCollector;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmomaticPE extends Trimmomatic
{
	private Logger logger;

	public TrimmomaticPE(Logger logger)
	{	
		this.logger=logger;
	}

	public void processPipeline(FastqParser rawParser1, FastqParser rawParser2, File output1P,File output1U, File output2P, File output2U, 
			Trimmer trimmers[], File trimLog, File statsSummary, PairingValidator pairingValidator, Boolean compressBlock, Integer compressLevel, int threads) throws Exception
	{
		boolean useParserWorkers = threads > 1;
		boolean useSerializerWorkers = threads > 1;
		boolean useParallelCompressor = compressBlock!=null? compressBlock: threads > 1;
	
		boolean useStatsWorker = threads > 1;
		boolean useLogWorker = threads > 1;
	
		ExceptionHolder exceptionHolder = new ExceptionHolder();
		
		Parser parser1 = Parser.makeParser(useParserWorkers, threads, rawParser1, exceptionHolder);
		Parser parser2 = Parser.makeParser(useParserWorkers, threads, rawParser2, exceptionHolder);

		Pipeline pipeline = Pipeline.makePipeline(threads, exceptionHolder);		
		
		Serializer serializer1P = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor, compressLevel, threads, output1P, exceptionHolder);
		Serializer serializer1U = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor, compressLevel, threads, output1U, exceptionHolder);
		Serializer serializer2P = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor, compressLevel, threads, output2P, exceptionHolder);
		Serializer serializer2U = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor, compressLevel, threads, output2U, exceptionHolder);

		List<Serializer> serializers = new ArrayList<Serializer>();
		serializers.add(serializer1P);
		serializers.add(serializer1U);
		serializers.add(serializer2P);
		serializers.add(serializer2U);
		
		TrimStatsCollector statsCollector = TrimStatsCollector.makeTrimStatsCollector(useStatsWorker, threads, exceptionHolder);
		TrimLogCollector logCollector = TrimLogCollector.makeTrimLogCollector(useLogWorker, threads, trimLog, exceptionHolder);
		
		boolean done = false;

		List<FastqRecord> recs1 = null;
		List<FastqRecord> recs2 = null;

		while (!done)
			{
			boolean done1 = false, done2 = false;
				
			if (!done1)
				{
				recs1 = null;
				while (recs1 == null)
					recs1 = parser1.poll();
				
				done1 = recs1.size()==0;
				}
			if (!done2)
				{
				recs2 = null;
				while (recs2 == null)
					recs2 = parser2.poll();

				done2 = recs2.size()==0;
				}

			done = done1 && done2;
				
			if(pairingValidator!=null)
				pairingValidator.validatePairs(recs1, recs2);
				
			BlockOfRecords bor = new BlockOfRecords(recs1, recs2);
			BlockOfWork work = new BlockOfWork(logger, trimmers, bor, done, true, trimLog != null, serializers, exceptionHolder);

			List<SerializedBlock> buffers = work.getBlocks();
				
			serializer1P.queueForWrite(buffers.get(0), exceptionHolder);
			serializer1U.queueForWrite(buffers.get(1), exceptionHolder);
			serializer2P.queueForWrite(buffers.get(2), exceptionHolder);
			serializer2U.queueForWrite(buffers.get(3), exceptionHolder);
								
			Future<BlockOfRecords> future = pipeline.submit(work);

			serializer1P.pollWritable();
			serializer1U.pollWritable();
			serializer2P.pollWritable();
			serializer2U.pollWritable();
				
			statsCollector.put(future);

			if (logCollector != null)
				logCollector.put(future);
			}
			
		parser1.close();
		parser2.close();
			
		pipeline.close();
		
		serializer1P.close();
		serializer1U.close();
		serializer2P.close();
		serializer2U.close();

		if (logCollector != null)
			logCollector.close();

		statsCollector.close();
			
		logger.infoln(statsCollector.getStats().processStatsPE(statsSummary));
	}

	public void process(File input1, File input2, File output1P, File output1U, File output2P, File output2U,
			Trimmer trimmers[], int phredOffset, File trimLog, File statsSummary, boolean validatePairing, Boolean compressBlock, Integer compressLevel, int threads) throws Exception
	{
		FastqParser parser1 = new FastqParser(phredOffset);
		parser1.open(input1);

		FastqParser parser2 = new FastqParser(phredOffset);
		parser2.open(input2);

		if(phredOffset==0)
			{
			int phred1=parser1.determinePhredOffset();
			int phred2=parser2.determinePhredOffset();
			
			if(phred1==phred2 && phred1!=0)
				{
				logger.infoln("Quality encoding detected as phred"+phred1);
				parser1.setPhredOffset(phred1);
				parser2.setPhredOffset(phred1);
				}
			else
				{
				logger.errorln("Error: Unable to detect quality encoding");
				System.exit(1);
				}
			}
		
		PairingValidator pairingValidator=null;
		
		if(validatePairing)
			pairingValidator=new PairingValidator(logger);
		
		processPipeline(parser1, parser2, output1P, output1U, output2P, output2U, trimmers,
				trimLog, statsSummary, pairingValidator, compressBlock, compressLevel, threads);

	}
	
	private static int getFileExtensionIndex(String str)
	{
		String extensions[]={".fq",".fastq",".txt",".gz",".bz2",".zip"};
	
		String tmp=str;
		boolean done=false;
		
		while(!done)
			{
			done=true;
			for(String ext: extensions)
				{		
				if(tmp.endsWith(ext))
					{
					tmp=tmp.substring(0,tmp.length()-ext.length());
					done=false;
					}
				}
			}
	
		return tmp.length();
	}

	private static String replaceLast(String str, String out, String in)
	{
		int idx1=str.lastIndexOf(out);
		if(idx1==-1)
			return null;
		
		int idx2=idx1+out.length();
		
		return str.substring(0,idx1)+in+str.substring(idx2);
	}
	
	
	
	private static File[] calculateTemplatedInput(String baseStr)
	{
		String translation[][]={{"_R1_","_R2_"},{"_f","_r"},{".f",".r"},{"_1","_2"},{".1",".2"}};
	
		File fileBase=new File(baseStr);
		File baseDir=fileBase.getParentFile();
		
		String baseName=fileBase.getName();
		int extSplit=getFileExtensionIndex(baseName);
		
		String core=baseName.substring(0,extSplit);
		String exts=baseName.substring(extSplit);
		
		for(String pair[]: translation)
			{
			String tmp=replaceLast(core, pair[0], pair[1]);
			if(tmp!=null)
				return new File[] {fileBase, new File(baseDir, tmp+exts)};
			}
		
		return null;
	}

	
	private static File[] calculateTemplatedOutput(String baseStr)
	{
		File fileBase=new File(baseStr);
		File baseDir=fileBase.getParentFile();
		
		String baseName=fileBase.getName();
		int extSplit=getFileExtensionIndex(baseName);
		
		String core=baseName.substring(0,extSplit);
		String exts=baseName.substring(extSplit);
		
		return new File[] {new File(baseDir,core+"_1P"+exts),new File(baseDir,core+"_1U"+exts),new File(baseDir,core+"_2P"+exts),new File(baseDir,core+"_2U"+exts)};		
	}

	

	public static boolean run(String[] args) throws Exception
	{
		int argIndex = 0;
		int phredOffset = 0;
		int threads = 0;

		boolean badOption = false;

		File trimLog = null;
		File statsSummary = null;
		
		boolean validatePairs = false;		
		boolean quiet=false;
		boolean showVersion=false;		
		
		Boolean compressBlock=null;
		Integer compressLevel=null;
		
		String templateInput=null;
		String templateOutput=null;
	

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
				else if (arg.equals("-basein"))
					{
					if (argIndex < args.length)
						templateInput = args[argIndex++];  
					else
						badOption = true;
					}
				else if (arg.equals("-baseout"))
					{
					if (argIndex < args.length)
						templateOutput = args[argIndex++];
					else
						badOption = true;
					}
				else if (arg.equals("-validatePairs"))
					validatePairs=true;
				else if (arg.equals("-compressLevel"))
					{
					if (argIndex < args.length)
						{
						compressLevel = Integer.parseInt(args[argIndex++]);
						if (compressLevel < 1 || compressLevel > 9)
							{
							System.err.println("compressLevel '"+compressLevel+"' should be between 1 and 9 inclusive");
							badOption = true;
							}
						}
					else
						badOption = true;
					}
				else if (arg.equals("-compressStream"))
					compressBlock=false;
				else if (arg.equals("-compressBlock"))
					compressBlock=true;	
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
		
		
		int additionalArgs=1+(templateInput==null?2:0)+(templateOutput==null?4:0);
		
		if ((nonOptionArgs.size() < additionalArgs) || badOption)
			return showVersion;
		
		Logger logger=new Logger(true,true,!quiet);
		
		
		logger.infoln("TrimmomaticPE: Started with arguments:");
		for (String arg : args)
			logger.info(" " + arg);
		logger.infoln();
		
		if(threads==0)
			{
			threads=calcAutoThreadCount();
			if(threads>1)
				logger.infoln("Multiple cores found: Using "+threads+" threads");
			}

		
		Iterator<String> nonOptionArgsIter=nonOptionArgs.iterator();

		File inputs[],outputs[];
		
		if(templateInput!=null)
			{
			inputs=calculateTemplatedInput(templateInput);
			if(inputs==null)
				{
				logger.errorln("Unable to determine input files from: "+templateInput);
				System.exit(1);
				}
			
			logger.infoln("Using templated Input files: "+inputs[0]+" "+inputs[1]);
			}
		else
			{
			inputs=new File[2];
			inputs[0]=new File(nonOptionArgsIter.next());
			inputs[1]=new File(nonOptionArgsIter.next());
			}
		
		if(templateOutput!=null)
			{
			outputs=calculateTemplatedOutput(templateOutput);
			if(outputs==null)
				{
				System.err.println("Unable to determine output files from: "+templateInput);
				System.exit(1);
				}
			
			logger.infoln("Using templated Output files: "+outputs[0]+" "+outputs[1]+" "+outputs[2]+" "+outputs[3]);
			}
		else
			{
			outputs=new File[4];
			outputs[0]=new File(nonOptionArgsIter.next());
			outputs[1]=new File(nonOptionArgsIter.next());
			outputs[2]=new File(nonOptionArgsIter.next());
			outputs[3]=new File(nonOptionArgsIter.next());
			}

		Trimmer trimmers[]=createTrimmers(logger, nonOptionArgsIter);
				
		TrimmomaticPE tm = new TrimmomaticPE(logger);
		tm.process(inputs[0], inputs[1], outputs[0], outputs[1], outputs[2], outputs[3], 
				trimmers, phredOffset, trimLog, statsSummary, validatePairs, compressBlock, compressLevel, threads);

		logger.infoln("TrimmomaticPE: Completed successfully");
		return true;
	}

	public static void main(String[] args) throws Exception
	{
		if (!run(args))
			{
			System.err.println("Usage: [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-compressLevel <lvl>] [-compressStream|-compressBlock] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.exit(1);
			}
	}

}
