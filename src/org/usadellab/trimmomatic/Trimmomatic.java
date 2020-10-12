package org.usadellab.trimmomatic;

import java.io.IOException;
import java.util.Arrays;

public class Trimmomatic
{

	private static final int MAX_AUTO_THREADS=16;

	public static int calcAutoThreadCount()
	{
		int cpus=Runtime.getRuntime().availableProcessors();
		
		if(cpus>MAX_AUTO_THREADS)
			return MAX_AUTO_THREADS;
		
		return cpus;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException
	{
		boolean showUsage=true;
	
		if(args.length>1)
			{
			String mode=args[0];
			String restOfArgs[]=Arrays.copyOfRange(args, 1, args.length);

			if(mode.equals("PE"))
				{
				if(TrimmomaticPE.run(restOfArgs))
					showUsage=false;
				}
			else if(mode.equals("SE"))
				{
				if(TrimmomaticSE.run(restOfArgs))
					showUsage=false;
				}
			}
		
		if(showUsage)
			{
			System.err.println("Usage: ");
			System.err.println("       PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.err.println("   or: ");
			System.err.println("       SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] <inputFile> <outputFile> <trimmer1>...");
			System.exit(1);
			}
	}

}
