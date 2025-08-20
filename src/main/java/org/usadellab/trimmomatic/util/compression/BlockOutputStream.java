package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.OutputStream;

public class BlockOutputStream {
	OutputStream stream;

	public BlockOutputStream(OutputStream stream) {
		this.stream = stream;
	}

	public void writeBytes(byte[] data) throws IOException {
		stream.write(data);
	}

	public void writeBlock(BlockData block) throws IOException {
		stream.write(block.getData());
	}

	public void close() throws IOException {
		stream.close();
	}

}
