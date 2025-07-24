package org.usadellab.trimmomatic.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;

import org.itadaki.bzip2.BZip2OutputStream;
import org.usadellab.trimmomatic.util.compression.CompressionFormat;

public class FastqSerializer {

	private BufferedWriter writer;
	private File inputFile;

	public FastqSerializer()
	{

	}

	public void open(OutputStream contentStream) throws IOException
	{
		writer = new BufferedWriter(new OutputStreamWriter(contentStream), 32768);
	}
	
	public void open(File file) throws IOException
	{
		this.inputFile = file;

		OutputStream rawStream = new FileOutputStream(file);		
		OutputStream contentStream = CompressionFormat.wrapStreamForSerializing(rawStream, file.getName(), null);
		
		open(contentStream);
	}

	public void close() throws IOException
	{
		writer.close();
	}

	public void writeRecord(FastqRecord record) throws IOException
	{	
		StringBuilder sb=new StringBuilder(record.getRecordLength());
			
	    sb.append('@');
	    sb.append(record.getName());
	    sb.append('\n');
	    sb.append(record.getSequence());
	    sb.append("\n+");
	    sb.append(record.getComment());
	    sb.append('\n');
	    sb.append(record.getQuality());
	    sb.append('\n');

	    writer.write(sb.toString());
	}

	public File getInputFile()
	{
		return inputFile;
	}

}
