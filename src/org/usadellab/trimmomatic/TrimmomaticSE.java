package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Future;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
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

public class TrimmomaticSE extends Trimmomatic {
	private Logger logger;

	public TrimmomaticSE(Logger logger) {
		this.logger = logger;
	}


	public void processPipeline(FastqParser rawParser, File output, Trimmer trimmers[],
								File trimLog, File statsSummary, Boolean compressBlock,
								Integer compressLevel, int threads) throws Exception {
		boolean useParserWorker = threads > 1;
		boolean useSerializerWorker = threads > 1;
		boolean useParallelCompressor = compressBlock != null ? compressBlock : threads > 1;
		
		boolean useStatsWorker = threads > 1;
		boolean useLogWorker = threads > 1;
		
		ExceptionHolder exceptionHolder = new ExceptionHolder();
		
		Parser parser = Parser.makeParser(useParserWorker, threads, rawParser, exceptionHolder);
		Pipeline pipeline = Pipeline.makePipeline(threads, exceptionHolder);		
		Serializer serializer = Serializer.makeSerializer(logger, useSerializerWorker, useParallelCompressor,
															compressLevel, threads, output, exceptionHolder);

		List<Serializer> serializers = new ArrayList<Serializer>();
		serializers.add(serializer);
		
		TrimStatsCollector statsCollector = TrimStatsCollector.makeTrimStatsCollector(useStatsWorker, threads, exceptionHolder);
		TrimLogCollector logCollector = TrimLogCollector.makeTrimLogCollector(useLogWorker, threads, trimLog, exceptionHolder);

		boolean done = false;

		List<FastqRecord> recs1 = null;

		while (!done) {
			recs1 = null;
			while (recs1 == null) recs1 = parser.poll();

			done = recs1.isEmpty();

			BlockOfRecords bor = new BlockOfRecords(recs1, null);
			BlockOfWork work = new BlockOfWork(logger, trimmers, bor, done, false, trimLog != null, serializers, exceptionHolder);

			List<SerializedBlock> buffers = work.getBlocks();
			
			serializer.queueForWrite(buffers.getFirst(), exceptionHolder);
			
			Future<BlockOfRecords> future = pipeline.submit(work);

			serializer.pollWritable();
			
			statsCollector.put(future);

			if (logCollector != null) {
				logger.infoln("queue for log");
				logCollector.put(future);
			}
		}

		parser.close();	
		pipeline.close();	
		serializer.close();
			
		if (logCollector != null) logCollector.close();

		statsCollector.close();
		logger.infoln(statsCollector.getStats().processStatsSE(statsSummary));
	}

	public void process(File input, File output, Trimmer trimmers[], int phredOffset, File trimLog,
						File statsSummary, Boolean compressBlock, Integer compressLevel, int threads) throws Exception {
		FastqParser parser = new FastqParser(phredOffset);
		parser.open(input);

		if(phredOffset == 0) {
			int phred = parser.determinePhredOffset();
			if(phred != 0) {
				logger.infoln("Quality encoding detected as phred" + phred);
				parser.setPhredOffset(phred);
			} else {
				logger.errorln("Error: Unable to detect quality encoding");
				System.exit(1);
			}
		}
		
		processPipeline(parser, output, trimmers, trimLog, statsSummary, compressBlock, compressLevel, threads);
	}

	public static boolean run(String[] args) throws Exception {
		int argIndex = 0;
		int phredOffset = 0;
		int threads = 0;

		boolean badOption = false;

		File trimLog = null;
		File statsSummary = null;
		
		boolean quiet = false;
		boolean showVersion = false;
		
		Boolean compressBlock = null;
		Integer compressLevel = null;
		
		List<String> nonOptionArgs = new ArrayList<String>();
		
		while (argIndex < args.length) {
			String arg = args[argIndex++];
			
			if(arg.startsWith("-")) {
				if (arg.equals("-phred33")) phredOffset = 33;
				else if (arg.equals("-phred64")) phredOffset = 64;
				else if (arg.equals("-threads")) {
					if (argIndex < args.length) threads = Integer.parseInt(args[argIndex++]);
					else badOption = true;
				}
				else if (arg.equals("-trimlog")) {
					if (argIndex < args.length) trimLog = new File(args[argIndex++]);
					else badOption = true;
				}
				else if (arg.equals("-summary")) {
					if (argIndex < args.length) statsSummary = new File(args[argIndex++]);
					else badOption = true;
				}
				else if (arg.equals("-compressLevel")) {
					if (argIndex < args.length) {
						compressLevel = Integer.parseInt(args[argIndex++]);
						if (compressLevel < 1 || compressLevel > 9) {
							System.err.println("compressLevel '" + compressLevel + "' should be between 1 and 9 inclusive");
							badOption = true;
						}
					} else badOption = true;
				}
				else if (arg.equals("-compressStream")) compressBlock = false;
				else if (arg.equals("-compressBlock")) compressBlock = true;
				else if (arg.equals("-quiet")) quiet = true;
				else if (arg.equals("-version")) showVersion = true;
				else {
					System.err.println("Unknown option " + arg);
					badOption = true;
				}
			} else nonOptionArgs.add(arg);
		}

		if (showVersion) Trimmomatic.showVersion();
		
		if ((nonOptionArgs.size() < 3) || badOption) return showVersion;

		Logger logger = new Logger(true, true, !quiet);
		
		logger.infoln("TrimmomaticSE: Started with arguments:");
		for (String arg : args) logger.info(" " + arg);
		logger.infoln();
		
		if(threads == 0) {
			threads = calcAutoThreadCount();
			logger.infoln("Automatically using " + threads + " threads");
		}
		
		Iterator<String> nonOptionArgsIter = nonOptionArgs.iterator();
		
		File input = new File(nonOptionArgsIter.next());
		File output = new File(nonOptionArgsIter.next());

		Trimmer[] trimmers = createTrimmers(logger, nonOptionArgsIter);

		TrimmomaticSE tm = new TrimmomaticSE(logger);
		tm.process(input, output, trimmers, phredOffset, trimLog, statsSummary, compressBlock, compressLevel, threads);

		logger.infoln("TrimmomaticSE: Completed successfully");
		return true;
	}

	public static void main(String[] args) throws Exception {
		if (!run(args)) {
			System.err.println("Usage: TrimmomaticSE [-version] [-threads <threads>] [-phred33|-phred64] " +
					"[-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-compressLevel <lvl>] " +
					"[-compressStream|-compressBlock] <inputFile> <outputFile> <trimmer1>...");

			System.exit(1);
		}
	}

}
