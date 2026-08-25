package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.trim.TrimmerFactory;
import org.usadellab.trimmomatic.util.Logger;

public class Trimmomatic {
	private static final int MAX_AUTO_THREADS_THRESHOLD = 8;
	private static final int MAX_AUTO_THREADS_ALLOC = 4;

	static void showVersion() {
		try {
			InputStream is = ClassLoader.getSystemResourceAsStream("version.properties");

			Properties props = new Properties();
			props.load(is);

			String version = props.getProperty("version");
			String hash    = props.getProperty("build.hash", "");
			String dirty   = props.getProperty("build.dirty", "false");

			// The placeholder is only substituted when the git-commit-id plugin runs
			// successfully.  If it didn't (source tarball, no git history, etc.) the
			// raw Maven placeholder string is left in the file, skip it in that case.
			StringBuilder sb = new StringBuilder(version);
			if (!hash.isEmpty() && !hash.startsWith("${")) {
				sb.append('+').append(hash);
				if ("true".equalsIgnoreCase(dirty)) {
					sb.append("-dirty");
				}
			}

			System.out.println(sb);
		} catch (Exception e) {
			throw new RuntimeException("Unable to determine version", e);
		}
	}

	static int calcAutoThreadCount() {
		int cpus = Runtime.getRuntime().availableProcessors();

		if (cpus <= MAX_AUTO_THREADS_THRESHOLD) {
			if (cpus < MAX_AUTO_THREADS_ALLOC)
				return cpus;

			return MAX_AUTO_THREADS_ALLOC;
		}

		return 1;
	}

	static Trimmer[] createTrimmers(Logger logger, Iterator<String> nonOptionArgsIter) throws IOException {
		TrimmerFactory fac = new TrimmerFactory(logger);

		List<Trimmer> trimmerList = new ArrayList<Trimmer>();
		while (nonOptionArgsIter.hasNext())
			trimmerList.add(fac.makeTrimmer(nonOptionArgsIter.next()));

		Trimmer trimmers[] = trimmerList.toArray(new Trimmer[0]);

		return trimmers;
	}

	/**
	 * Parses a single whitespace-separated step-list string (e.g. the value of
	 * -pe1steps/-pe2steps) into a Trimmer array. An empty/blank string yields a
	 * zero-length array (pure passthrough for that mate), matching how an empty
	 * trailing step list behaves.
	 */
	static Trimmer[] createTrimmersFromString(Logger logger, String stepsStr) throws IOException {
		if (stepsStr == null || stepsStr.trim().isEmpty())
			return new Trimmer[0];

		List<String> tokens = Arrays.asList(stepsStr.trim().split("\\s+"));
		return createTrimmers(logger, tokens.iterator());
	}

	static String trimExtension(String filename) {
		String extensions[] = { ".fq", ".fastq", ".txt", ".gz", ".bz2", ".zip" };

		String tmp = filename;
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

		return tmp;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		boolean showUsage = true;

		if (args.length > 0) {
			String mode = args[0];
			String restOfArgs[] = Arrays.copyOfRange(args, 1, args.length);

			if (mode.equals("PE")) {
				if (TrimmomaticPE.run(restOfArgs))
					showUsage = false;
			} else if (mode.equals("SE")) {
				if (TrimmomaticSE.run(restOfArgs))
					showUsage = false;
			} else if (mode.equals("-version")) {
				showVersion();
				showUsage = false;
			} else if (args.length == 1) {
				File input = new File(args[0]);
				String inputPath = input.getAbsolutePath();
				String output = trimExtension(inputPath) + ".trimmed.fq.gz";

				List<String> seArgs = new ArrayList<String>();
				seArgs.add("-threads");
				seArgs.add(String.valueOf(Runtime.getRuntime().availableProcessors()));
				seArgs.add(inputPath);
				seArgs.add(output);
				seArgs.add("ILLUMINACLIP:TruSeq3-SE-GGGGG.fa:2:30:10");
				seArgs.add("SLIDINGWINDOW:4:20");
				seArgs.add("MINLEN:36");

				if (TrimmomaticSE.run(seArgs.toArray(new String[0])))
					showUsage = false;
			} else if (args.length == 2) {
				File input1 = new File(args[0]);
				File input2 = new File(args[1]);
				String input1Path = input1.getAbsolutePath();
				String input2Path = input2.getAbsolutePath();
				String base1 = trimExtension(input1Path);
				String base2 = trimExtension(input2Path);

				List<String> peArgs = new ArrayList<String>();
				peArgs.add("-threads");
				peArgs.add(String.valueOf(Runtime.getRuntime().availableProcessors()));
				peArgs.add(input1Path);
				peArgs.add(input2Path);
				peArgs.add(base1 + ".trimmed.paired.fq.gz");
				peArgs.add(base1 + ".trimmed.unpaired.fq.gz");
				peArgs.add(base2 + ".trimmed.paired.fq.gz");
				peArgs.add(base2 + ".trimmed.unpaired.fq.gz");
				peArgs.add("ILLUMINACLIP:TruSeq3-PE-2-GGGGG.fa:2:30:10");
				peArgs.add("SLIDINGWINDOW:4:20");
				peArgs.add("MINLEN:36");

				if (TrimmomaticPE.run(peArgs.toArray(new String[0])))
					showUsage = false;
			}
		}

		if (showUsage) {
			System.err.println("Usage: ");
			System.err.println(
					"       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.err.println("   or: ");
			System.err.println(
					"       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...");
			System.err.println("   or: ");
			System.err.println("       -version");
			System.err.println("   or: ");
			System.err.println("       <inputFile1> [inputFile2] (for automatic parameter selection)");
			System.exit(1);
		}
	}

}
