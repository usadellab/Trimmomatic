package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;

public class Pairomatic
{

	/**
	 * Pairomatic: The FASTQ pair/unpairer
	 */

	public Pairomatic()
	{

	}

	private Set<String> getFastqNames(File file, Character delimiter) throws IOException
	{
		Set<String> names = new LinkedHashSet<String>();

		FastqParser parser = new FastqParser(0);
		parser.parse(file);

		while (parser.hasNext())
			{
			FastqRecord rec = parser.next();

			String name = rec.getName();

			if (delimiter != null)
				{
				int index = name.lastIndexOf(delimiter);

				if (index == -1)
					throw new RuntimeException("Error: Failed to find expected delimiter '" + delimiter
							+ "' in record named '" + name + "'");

				name = name.substring(0, index);
				}
			
			if(names.contains(name))
				throw new RuntimeException("Error: Found "+name+" more than once in file - check delimiter is correct '"+delimiter+"'");
			
			names.add(name);
			}

		return names;
	}

	
	private boolean equalOrdering(Set<String> set1, Set<String> set2)
	{
		if(set1.size()!=set2.size())
			return false;
		
		Iterator<String> iter1=set1.iterator();
		Iterator<String> iter2=set2.iterator();
		
		while(iter1.hasNext() && iter2.hasNext())
			{
			String str1=iter1.next();
			String str2=iter2.next();
			
			if(!str1.equals(str2))
				return false;
			}
		
		if(iter1.hasNext())
			return false;
		
		if(iter2.hasNext())
			return false;
		
		return true;
	}
	
	
	private void splitFastq(File input, File match, File unmatch, Set<String> toKeep, Character delimiter)
			throws IOException
	{
		FastqParser parser = new FastqParser(0);
		parser.parse(input);

		FastqSerializer matchSerializer=new FastqSerializer();
		matchSerializer.open(match);
		
		FastqSerializer unmatchSerializer=new FastqSerializer();
		unmatchSerializer.open(unmatch);
		
		while (parser.hasNext())
			{
			FastqRecord rec = parser.next();

			String name = rec.getName();

			if (delimiter != null)
				{
				int index = name.indexOf(delimiter);

				if (index == -1)
					throw new RuntimeException("Failed to find expected delimiter '" + delimiter
							+ "' in record named '" + name + "'");

				name = name.substring(0, index);
				}
			
			if(toKeep.contains(name))
				matchSerializer.writeRecord(rec);
			else
				unmatchSerializer.writeRecord(rec);
			}
		
		matchSerializer.close();
		unmatchSerializer.close();

	}

	public void process(File input1, File input2, File output1P, File output1U, File output2P, File output2U,
			Character delimiter) throws IOException
	{
		Set<String> names1 = getFastqNames(input1, delimiter);
		System.err.println("First input file contains " + names1.size() + " records");

		Set<String> names2 = getFastqNames(input2, delimiter);
		System.err.println("Second input file contains " + names2.size() + " records");

		names1.retainAll(names2);
		System.err.println("Files shared " + names1.size() + " records");
		
		names2.retainAll(names1);
		if(!equalOrdering(names1,names2))
			{
			System.err.println("Error: Common records are not in identical order, cowardly refusing to do anything");
			return;
			}
		
		System.err.println("Splitting first file");
		splitFastq(input1, output1P, output1U, names1, delimiter);
		
		System.err.println("Splitting second file");
		splitFastq(input2, output2P, output2U, names1, delimiter);
		
		System.err.println("All done");
	}

	public static void main(String[] args) throws IOException
	{
		int argIndex = 0;
		Character delim=null;
		
		boolean badOption = false;

		while (argIndex < args.length && args[argIndex].startsWith("-"))
			{
			String arg = args[argIndex++];
			if(arg.equals("-delim"))
				{
				String delimStr=args[argIndex++];
				
				if(delimStr.length()!=1)
					System.err.println("Delimiter must be exactly one character, got '"+delimStr+"'");
				else
					delim=delimStr.charAt(0);
				}
			else
				{
				System.err.println("Unknown option " + arg);
				badOption = true;
				}
			}

		if (args.length - argIndex < 6 || badOption)
			{
			System.err
					.println("Usage: Pairomatic [-delim delimChar] <inputFile1> <inputFile2> <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>");
			System.exit(1);
			}		
		
		File input1 = new File(args[argIndex++]);
		File input2 = new File(args[argIndex++]);

		File output1P = new File(args[argIndex++]);
		File output1U = new File(args[argIndex++]);

		File output2P = new File(args[argIndex++]);
		File output2U = new File(args[argIndex++]);

		Pairomatic pm = new Pairomatic();
		pm.process(input1, input2, output1P, output1U, output2P, output2U, delim);
	
	}
}
