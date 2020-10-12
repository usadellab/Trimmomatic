package org.usadellab.trimmomatic.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;

import org.itadaki.bzip2.BZip2OutputStream;

public class FastqSerializer {

	private BufferedWriter stream;
	private File inputFile;

	public FastqSerializer()
	{

	}

	public void open(File file) throws IOException
	{
		String name = file.getName();
		this.inputFile = file;

		OutputStream gStream = new FileOutputStream(file);

		if (name.endsWith(".gz"))
			{
			gStream = new GZIPOutputStream(gStream);
			}
		else if (name.endsWith(".bz2"))
			{
			gStream = new BZip2OutputStream(gStream);
			}

		// stream=new OutputStreamWriter(new BufferedOutputStream(gStream));

		stream = new BufferedWriter(new OutputStreamWriter(gStream), 32768);
	}

	public void close() throws IOException
	{
		stream.close();
	}

	public void writeRecord(FastqRecord record) throws IOException
	{
		StringBuilder sb=new StringBuilder(500);
	       
	    sb.append('@');
	    sb.append(record.getName());
	    sb.append('\n');
	    sb.append(record.getSequence());
	    sb.append("\n+");
	    sb.append(record.getComment());
	    sb.append('\n');
	    sb.append(record.getQuality());
	    sb.append('\n');

	    stream.write(sb.toString());
	}

	public File getInputFile()
	{
		return inputFile;
	}

	public void setInputFile(File file)
	{
		this.inputFile = file;
	}

}
