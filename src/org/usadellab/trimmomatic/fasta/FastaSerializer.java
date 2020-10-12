package org.usadellab.trimmomatic.fasta;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

public class FastaSerializer {

	public static final int BASES_PER_LINE=60;
	
	private PrintStream stream;
	
	public FastaSerializer()
	{
		
	}
	
	public void open(File file) throws IOException
	{
		stream=new PrintStream(file);
	}

	public void open(PrintStream stream) throws IOException
	{
		this.stream=new PrintStream(stream);
	}
	
	public void close() throws IOException
	{
		stream.close();
	}
	
	public synchronized void writeRecord(FastaRecord record) 
	{
		stream.println(">"+record.getName());
		
		String seq=record.getSequence();
		
		int len=seq.length();
		int pos=0;
		while(pos+BASES_PER_LINE <= len)
			{
			stream.println(seq.substring(pos,pos+BASES_PER_LINE));
			pos+=BASES_PER_LINE;
			}

		String remainder=seq.substring(pos);
		if(remainder.length()>0)
			stream.println(remainder);
	}
	
}
