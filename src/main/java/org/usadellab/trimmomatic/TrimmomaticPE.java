package org.usadellab.trimmomatic;

import java.io.File;
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
import org.usadellab.trimmomatic.threading.parser.InterleavedParserPair;
import org.usadellab.trimmomatic.threading.parser.Parser;
import org.usadellab.trimmomatic.threading.pipeline.Pipeline;
import org.usadellab.trimmomatic.threading.serializer.SerializedBlock;
import org.usadellab.trimmomatic.threading.serializer.Serializer;
import org.usadellab.trimmomatic.threading.trimlog.TrimLogCollector;
import org.usadellab.trimmomatic.threading.trimstats.TrimStatsCollector;
import org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmomaticPE extends Trimmomatic {
	private Logger logger;

	public TrimmomaticPE(Logger logger) {
		this.logger = logger;
	}

	public void processPipeline(FastqParser rawParser1, FastqParser rawParser2, boolean interleaved, File output1P,
			File output1U, File output2P, File output2U, Trimmer trimmers[], File trimLog, File statsSummary,
			PairingValidator pairingValidator, Boolean compressBlock, Integer compressLevel, int threads, boolean verbose,
			int technicalRead) throws Exception {
		processPipeline(rawParser1, rawParser2, interleaved, output1P, output1U, output2P, output2U, trimmers, null,
				null, trimLog, statsSummary, pairingValidator, compressBlock, compressLevel, threads, verbose,
				technicalRead);
	}

	/**
	 * trimmers1/trimmers2 (both non-null together, or both null) select per-mate
	 * step routing: mate 1 runs only trimmers1, mate 2 runs only trimmers2, and
	 * either mate's steps returning null drops the whole pair. Mutually exclusive
	 * with technicalRead != 0, see BlockOfWork's constructor validation.
	 */
	public void processPipeline(FastqParser rawParser1, FastqParser rawParser2, boolean interleaved, File output1P,
			File output1U, File output2P, File output2U, Trimmer trimmers[], Trimmer trimmers1[], Trimmer trimmers2[],
			File trimLog, File statsSummary, PairingValidator pairingValidator, Boolean compressBlock,
			Integer compressLevel, int threads, boolean verbose, int technicalRead) throws Exception {
		boolean useParserWorkers = threads > 1;
		boolean useSerializerWorkers = threads > 1;
		boolean useParallelCompressor = compressBlock != null ? compressBlock : threads > 1;

		boolean useStatsWorker = threads > 1;
		boolean useLogWorker = threads > 1;

		ExceptionHolder exceptionHolder = new ExceptionHolder();

		final InterleavedParserPair interleavedPair;
		final Parser parser1, parser2;

		if (interleaved) {
			interleavedPair = new InterleavedParserPair(rawParser1, threads, exceptionHolder);
			parser1 = interleavedPair.getR1Parser();
			parser2 = interleavedPair.getR2Parser();
		} else {
			interleavedPair = null;
			parser1 = Parser.makeParser(useParserWorkers, threads, rawParser1, exceptionHolder);
			parser2 = Parser.makeParser(useParserWorkers, threads, rawParser2, exceptionHolder);
		}

		// Resources are declared in the reverse of the order they must be closed in
		// (parser1/parser2 first, statsCollector last), try-with-resources closes
		// bottom-to-top, so cleanup still happens if the loop below throws, or if a
		// later resource here fails to construct after an earlier one already opened
		// a file / started a background thread.
		// parser1/parser2 are no-op HalfParsers when interleaved (interleavedPair owns
		// the real source in that case), so closing them unconditionally is harmless.
		try (TrimStatsCollector statsCollector = TrimStatsCollector.makeTrimStatsCollector(useStatsWorker, threads,
				exceptionHolder);
				TrimLogCollector logCollector = TrimLogCollector.makeTrimLogCollector(useLogWorker, threads, trimLog,
						exceptionHolder);
				Serializer serializer2U = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor,
						compressLevel, threads, output2U, exceptionHolder);
				Serializer serializer2P = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor,
						compressLevel, threads, output2P, exceptionHolder);
				Serializer serializer1U = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor,
						compressLevel, threads, output1U, exceptionHolder);
				Serializer serializer1P = Serializer.makeSerializer(logger, useSerializerWorkers, useParallelCompressor,
						compressLevel, threads, output1P, exceptionHolder);
				Pipeline pipeline = Pipeline.makePipeline(threads, exceptionHolder);
				interleavedPair;
				parser2;
				parser1) {

			List<Serializer> serializers = new ArrayList<Serializer>();
			serializers.add(serializer1P);
			serializers.add(serializer1U);
			serializers.add(serializer2P);
			serializers.add(serializer2U);

			boolean done = false;

			List<FastqRecord> recs1 = null;
			List<FastqRecord> recs2 = null;

			while (!done) {
				boolean done1 = false, done2 = false;

				if (!done1) {
					recs1 = null;
					while (recs1 == null)
						recs1 = parser1.poll();

					done1 = recs1.size() == 0;
				}
				if (!done2) {
					recs2 = null;
					while (recs2 == null)
						recs2 = parser2.poll();

					done2 = recs2.size() == 0;
				}

				done = done1 && done2;

				if (pairingValidator != null)
					pairingValidator.validatePairs(recs1, recs2);

				BlockOfRecords bor = new BlockOfRecords(recs1, recs2);
				BlockOfWork work = new BlockOfWork(logger, trimmers, trimmers1, trimmers2, bor, done, true,
						technicalRead, trimLog != null, serializers, exceptionHolder);

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

			// statsCollector merges each block's stats asynchronously on its own
			// background thread (SelfThreadedTrimStatsCollector, used when threads > 1).
			// Reading .getStats() without first waiting for that thread to finish
			// draining races the merge, for a fast-to-read input the main thread can
			// win and log/write an undercount despite the FASTQ output itself being
			// complete and correct. close() waits for the thread to finish first.
			statsCollector.close();
			logger.infoln(statsCollector.getStats().processStatsPE(statsSummary));

			if (verbose) {
				for (Trimmer t : trimmers) {
					if (t instanceof IlluminaClippingTrimmer ict)
						ict.printStats(logger);
				}
			}
		}
	}

	public void process(File input1, File input2, boolean interleaved, File output1P, File output1U, File output2P,
			File output2U, Trimmer trimmers[], int phredOffset, File trimLog, File statsSummary, boolean validatePairing,
			Boolean compressBlock, Integer compressLevel, int threads, boolean verbose, int technicalRead)
			throws Exception {
		process(input1, input2, interleaved, output1P, output1U, output2P, output2U, trimmers, null, null,
				phredOffset, trimLog, statsSummary, validatePairing, compressBlock, compressLevel, threads, verbose,
				technicalRead);
	}

	/**
	 * trimmers1/trimmers2 (both non-null together, or both null) select per-mate
	 * step routing, see processPipeline's per-mate overload.
	 */
	public void process(File input1, File input2, boolean interleaved, File output1P, File output1U, File output2P,
			File output2U, Trimmer trimmers[], Trimmer trimmers1[], Trimmer trimmers2[], int phredOffset,
			File trimLog, File statsSummary, boolean validatePairing, Boolean compressBlock, Integer compressLevel,
			int threads, boolean verbose, int technicalRead) throws Exception {
		FastqParser parser1 = new FastqParser(phredOffset);
		FastqParser parser2 = interleaved ? null : new FastqParser(phredOffset);

		if (interleaved) {
			parser1.open(input1);
		} else {
			Exception[] openErrors = new Exception[2];
			Thread t1 = Thread.ofVirtual().start(() -> {
				try { parser1.open(input1); }
				catch (Exception e) { openErrors[0] = e; }
			});
			Thread t2 = Thread.ofVirtual().start(() -> {
				try { parser2.open(input2); }
				catch (Exception e) { openErrors[1] = e; }
			});
			t1.join();
			t2.join();
			if (openErrors[0] != null || openErrors[1] != null) {
				// One side may have opened successfully before the other failed,
				// close it so its file handle isn't leaked. Only close a parser whose
				// open() actually succeeded; close() on a never-opened FastqParser
				// would otherwise touch an unassigned reader field.
				if (openErrors[0] == null) {
					try { parser1.close(); } catch (Exception closeEx) { /* best-effort cleanup */ }
				}
				if (openErrors[1] == null) {
					try { parser2.close(); } catch (Exception closeEx) { /* best-effort cleanup */ }
				}
				throw openErrors[0] != null ? openErrors[0] : openErrors[1];
			}
		}

		if (phredOffset == 0) {
			if (interleaved) {
				int phred1 = parser1.determinePhredOffset();
				if (phred1 != 0) {
					logger.infoln("Quality encoding detected as phred" + phred1);
					parser1.setPhredOffset(phred1);
				} else {
					logger.errorln("Error: Unable to detect quality encoding");
					System.exit(1);
				}
			} else {
				int phred1 = parser1.determinePhredOffset();
				int phred2 = parser2.determinePhredOffset();

				if (phred1 == phred2 && phred1 != 0) {
					logger.infoln("Quality encoding detected as phred" + phred1);
					parser1.setPhredOffset(phred1);
					parser2.setPhredOffset(phred1);
				} else {
					logger.errorln("Error: Unable to detect quality encoding");
					System.exit(1);
				}
			}
		}

		PairingValidator pairingValidator = null;

		if (validatePairing && !interleaved)
			pairingValidator = new PairingValidator(logger);

		processPipeline(parser1, parser2, interleaved, output1P, output1U, output2P, output2U, trimmers, trimmers1,
				trimmers2, trimLog, statsSummary, pairingValidator, compressBlock, compressLevel, threads, verbose,
				technicalRead);

	}

	private static int getFileExtensionIndex(String str) {
		String extensions[] = { ".fq", ".fastq", ".txt", ".gz", ".bz2", ".zip" };

		String tmp = str;
		boolean done = false;

		while (!done) {
			done = true;
			for (String ext : extensions) {
				if (tmp.endsWith(ext)) {
					tmp = tmp.substring(0, tmp.length() - ext.length());
					done = false;
				}
			}
		}

		return tmp.length();
	}

	private static String replaceLast(String str, String out, String in) {
		int idx1 = str.lastIndexOf(out);
		if (idx1 == -1)
			return null;

		int idx2 = idx1 + out.length();

		return str.substring(0, idx1) + in + str.substring(idx2);
	}

	private static File[] calculateTemplatedInput(String baseStr) {
		String translation[][] = { { "_R1_", "_R2_" }, { "_f", "_r" }, { ".f", ".r" }, { "_1", "_2" }, { ".1", ".2" } };

		File fileBase = new File(baseStr);
		File baseDir = fileBase.getParentFile();

		String baseName = fileBase.getName();
		int extSplit = getFileExtensionIndex(baseName);

		String core = baseName.substring(0, extSplit);
		String exts = baseName.substring(extSplit);

		for (String pair[] : translation) {
			String tmp = replaceLast(core, pair[0], pair[1]);
			if (tmp != null)
				return new File[] { fileBase, new File(baseDir, tmp + exts) };
		}

		return null;
	}

	private static File[] calculateTemplatedOutput(String baseStr) {
		File fileBase = new File(baseStr);
		File baseDir = fileBase.getParentFile();

		String baseName = fileBase.getName();
		int extSplit = getFileExtensionIndex(baseName);

		String core = baseName.substring(0, extSplit);
		String exts = baseName.substring(extSplit);

		return new File[] { new File(baseDir, core + "_1P" + exts), new File(baseDir, core + "_1U" + exts),
				new File(baseDir, core + "_2P" + exts), new File(baseDir, core + "_2U" + exts) };
	}

	public static boolean run(String[] args) throws Exception {
		int argIndex = 0;
		int phredOffset = 0;
		int threads = 0;

		boolean badOption = false;

		File trimLog = null;
		File statsSummary = null;

		boolean validatePairs = false;
		boolean quiet = false;
		boolean showVersion = false;
		boolean verbose = false;
		boolean interleaved = false;
		boolean longread = false;
		int technicalRead = 0;
		String pe1StepsStr = null;
		String pe2StepsStr = null;

		Boolean compressBlock = null;
		Integer compressLevel = null;

		String templateInput = null;
		String templateOutput = null;

		List<String> nonOptionArgs = new ArrayList<String>();

		while (argIndex < args.length) {
			String arg = args[argIndex++];

			if (arg.startsWith("-")) {
				if (arg.equals("-phred33"))
					phredOffset = 33;
				else if (arg.equals("-phred64"))
					phredOffset = 64;
				else if (arg.equals("-threads"))
					threads = Integer.parseInt(args[argIndex++]);
				else if (arg.equals("-trimlog")) {
					if (argIndex < args.length)
						trimLog = new File(args[argIndex++]);
					else
						badOption = true;
				} else if (arg.equals("-summary")) {
					if (argIndex < args.length)
						statsSummary = new File(args[argIndex++]);
					else
						badOption = true;
				} else if (arg.equals("-basein")) {
					if (argIndex < args.length)
						templateInput = args[argIndex++];
					else
						badOption = true;
				} else if (arg.equals("-baseout")) {
					if (argIndex < args.length)
						templateOutput = args[argIndex++];
					else
						badOption = true;
				} else if (arg.equals("-validatePairs"))
					validatePairs = true;
				else if (arg.equals("-compressLevel")) {
					if (argIndex < args.length) {
						compressLevel = Integer.parseInt(args[argIndex++]);
						if (compressLevel < 1 || compressLevel > 9) {
							System.err.println(
									"compressLevel '" + compressLevel + "' should be between 1 and 9 inclusive");
							badOption = true;
						}
					} else
						badOption = true;
				} else if (arg.equals("-compressStream"))
					compressBlock = false;
				else if (arg.equals("-compressBlock"))
					compressBlock = true;
				else if (arg.equals("-quiet"))
					quiet = true;
				else if (arg.equals("-verbose"))
					verbose = true;
				else if (arg.equals("-version"))
					showVersion = true;
				else if (arg.equals("-interleaved"))
					interleaved = true;
				else if (arg.equals("-longread"))
					longread = true;
				else if (arg.equals("-technicalread")) {
					if (argIndex < args.length) {
						technicalRead = Integer.parseInt(args[argIndex++]);
						if (technicalRead != 1 && technicalRead != 2) {
							System.err.println("technicalread must be 1 or 2");
							badOption = true;
						}
					} else
						badOption = true;
				} else if (arg.equals("-pe1steps")) {
					if (argIndex < args.length)
						pe1StepsStr = args[argIndex++];
					else
						badOption = true;
				} else if (arg.equals("-pe2steps")) {
					if (argIndex < args.length)
						pe2StepsStr = args[argIndex++];
					else
						badOption = true;
				} else {
					System.err.println("Unknown option " + arg);
					badOption = true;
				}
			} else
				nonOptionArgs.add(arg);
		}

		boolean perMateMode = (pe1StepsStr != null) || (pe2StepsStr != null);

		if (perMateMode && (pe1StepsStr == null || pe2StepsStr == null)) {
			System.err.println("-pe1steps and -pe2steps must both be given together");
			badOption = true;
		}
		if (perMateMode && technicalRead != 0) {
			System.err.println("-pe1steps/-pe2steps cannot be combined with -technicalread");
			badOption = true;
		}

		if (showVersion)
			Trimmomatic.showVersion();

		int inputFiles = interleaved ? 1 : 2;
		// Normally at least one trailing trimmer step is required (the "+1"); in
		// per-mate mode the steps live in -pe1steps/-pe2steps instead, and the
		// trailing list must in fact be EMPTY (checked below), so that minimum
		// doesn't apply.
		int additionalArgs = (perMateMode ? 0 : 1) + (templateInput == null ? inputFiles : 0)
				+ (templateOutput == null ? 4 : 0);

		if ((nonOptionArgs.size() < additionalArgs) || badOption)
			return showVersion;

		Logger logger = new Logger(true, true, !quiet);

		logger.infoln("TrimmomaticPE: Started with arguments:");
		for (String arg : args)
			logger.info(" " + arg);
		logger.infoln();

		if (threads == 0) {
			threads = calcAutoThreadCount();
			if (threads > 1)
				logger.infoln("Multiple cores found: Using " + threads + " threads");
		}

		// -longread skips the expensive 10k-record phred pre-read and defaults to
		// Phred+33 (correct for all current long-read platforms: ONT, PacBio HiFi/CLR).
		if (longread && phredOffset == 0) {
			phredOffset = 33;
			logger.infoln("Long-read mode: assuming Phred+33 quality encoding");
		}

		Iterator<String> nonOptionArgsIter = nonOptionArgs.iterator();

		File inputs[], outputs[];

		if (interleaved) {
			// Interleaved mode: only one input file is needed.
			// -basein is treated as the direct path to the interleaved file, not as a template.
			String inputPath = (templateInput != null) ? templateInput : nonOptionArgsIter.next();
			inputs = new File[] { new File(inputPath), null };
			logger.infoln("Using interleaved input file: " + inputs[0]);
		} else if (templateInput != null) {
			inputs = calculateTemplatedInput(templateInput);
			if (inputs == null) {
				logger.errorln("Unable to determine input files from: " + templateInput);
				System.exit(1);
			}

			logger.infoln("Using templated Input files: " + inputs[0] + " " + inputs[1]);
		} else {
			inputs = new File[2];
			inputs[0] = new File(nonOptionArgsIter.next());
			inputs[1] = new File(nonOptionArgsIter.next());
		}

		if (templateOutput != null) {
			outputs = calculateTemplatedOutput(templateOutput);
			if (outputs == null) {
				System.err.println("Unable to determine output files from: " + templateInput);
				System.exit(1);
			}

			logger.infoln("Using templated Output files: " + outputs[0] + " " + outputs[1] + " " + outputs[2] + " "
					+ outputs[3]);
		} else {
			outputs = new File[4];
			outputs[0] = new File(nonOptionArgsIter.next());
			outputs[1] = new File(nonOptionArgsIter.next());
			outputs[2] = new File(nonOptionArgsIter.next());
			outputs[3] = new File(nonOptionArgsIter.next());
		}

		Trimmer trimmers[] = createTrimmers(logger, nonOptionArgsIter);

		if (perMateMode && trimmers.length > 0) {
			logger.errorln("Steps given both via -pe1steps/-pe2steps and as a trailing step list; "
					+ "put every step under one or the other, not both.");
			System.exit(1);
		}

		Trimmer trimmers1[] = null, trimmers2[] = null;
		if (perMateMode) {
			trimmers1 = createTrimmersFromString(logger, pe1StepsStr);
			trimmers2 = createTrimmersFromString(logger, pe2StepsStr);
		}

		TrimmomaticPE tm = new TrimmomaticPE(logger);
		tm.process(inputs[0], inputs[1], interleaved, outputs[0], outputs[1], outputs[2], outputs[3], trimmers,
				trimmers1, trimmers2, phredOffset, trimLog, statsSummary, validatePairs, compressBlock, compressLevel,
				threads, verbose, technicalRead);

		logger.infoln("TrimmomaticPE: Completed successfully");
		return true;
	}

	public static void main(String[] args) throws Exception {
		if (!run(args)) {
			System.err.println(
					"Usage: [-version] [-threads <threads>] [-phred33|-phred64] [-longread] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-verbose] [-validatePairs] [-interleaved] [-technicalread <1|2>] [-pe1steps <steps> -pe2steps <steps>] [-compressLevel <lvl>] [-compressStream|-compressBlock] [-basein <inputBase> | <inputFile1> [<inputFile2>]] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.exit(1);
		}
	}

}
