package org.usadellab.trimmomatic.util.compression;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;
//import org.usadellab.trimmomatic.threading.parser.Parser;

public class UncompressedBlockData implements BlockData {
	private byte[] data;

	public UncompressedBlockData(List<FastqRecord> recs) {
		ByteArrayOutputStream stream = new ByteArrayOutputStream();

		FastqSerializer ser = new FastqSerializer();

		try {
			ser.open(stream);

			for (FastqRecord rec : recs)
				ser.writeRecord(rec);

			ser.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		data = stream.toByteArray();
	}

	public byte[] getData() {
		return data;
	}
}
