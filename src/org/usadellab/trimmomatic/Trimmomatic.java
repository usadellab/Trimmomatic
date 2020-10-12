package org.usadellab.trimmomatic;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.trim.TrimmerFactory;
import org.usadellab.trimmomatic.util.Logger;

public class Trimmomatic
{

	private static final int MAX_AUTO_THREADS_THRESHOLD=8;
	private static final int MAX_AUTO_THREADS_ALLOC=4;

	static int calcAutoThreadCount()
	{
		int cpus=Runtime.getRuntime().availableProcessors();
		
		if(cpus<=MAX_AUTO_THREADS_THRESHOLD)
			{
			if(cpus<MAX_AUTO_THREADS_ALLOC)
				return cpus;
			
			return MAX_AUTO_THREADS_ALLOC;
			}
		
		return 1;
	}

	
	static Trimmer[] createTrimmers(Logger logger, Iterator<String> nonOptionArgsIter) throws IOException
	{
		TrimmerFactory fac = new TrimmerFactory(logger);
	
		List<Trimmer> trimmerList=new ArrayList<Trimmer>();
		while(nonOptionArgsIter.hasNext())
			trimmerList.add(fac.makeTrimmer(nonOptionArgsIter.next()));
	
		Trimmer trimmers[] = trimmerList.toArray(new Trimmer[0]);
		
		return trimmers;
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
			System.err.println("       PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.err.println("   or: ");
			System.err.println("       SE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...");
			System.exit(1);
			}
	}

}
