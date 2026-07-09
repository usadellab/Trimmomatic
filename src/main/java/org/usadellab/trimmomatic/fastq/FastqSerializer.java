package org.usadellab.trimmomatic.fastq;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import org.usadellab.trimmomatic.util.compression.CompressionFormat;

public class FastqSerializer implements Closeable {

	private BufferedWriter writer;
	private File inputFile;

	public FastqSerializer() {

	}

	public void open(OutputStream contentStream) throws IOException {
		writer = new BufferedWriter(new OutputStreamWriter(contentStream), 32768);
	}

	public void open(File file) throws IOException {
		this.inputFile = file;

		OutputStream rawStream = new FileOutputStream(file);
		OutputStream contentStream = CompressionFormat.wrapStreamForSerializing(rawStream, file.getName(), null);

		open(contentStream);
	}

	public void close() throws IOException {
		if (writer != null)
			writer.close();
	}

	public void writeRecord(FastqRecord record) throws IOException {
		writer.write('@');
		writer.write(record.getName());
		writer.write('\n');
		writer.write(record.getSequence());
		writer.write("\n+");
		writer.write(record.getComment());
		writer.write('\n');
		writer.write(record.getQuality());
		writer.write('\n');
	}

	public File getInputFile() {
		return inputFile;
	}

}
