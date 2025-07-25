package org.usadellab.trimmomatic.util.compression;

import java.io.ByteArrayOutputStream;
//import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.CRC32;
import java.util.zip.Deflater;
//import java.util.zip.DeflaterOutputStream;
//import java.util.zip.GZIPOutputStream;

import org.usadellab.trimmomatic.util.Logger;

public class GzipParallelCompressor implements ParallelCompressor {
	private static final byte GZIP_MAGIC_HIGH = (byte) 0x8b;
	private static final byte GZIP_MAGIC_LOW = (byte) 0x1f;

	@SuppressWarnings("unused")
	private Logger logger;
	private int compressLevel;

	private int length;
	private CRC32 crc;

	public GzipParallelCompressor(Logger logger, Integer compressLevel) {
		this.logger = logger;
		this.compressLevel = compressLevel == null ? Deflater.DEFAULT_COMPRESSION : compressLevel;

		crc = new CRC32();
		crc.reset();
	}

	@Override
	public BlockOutputStream wrapAndWriteHeader(OutputStream stream) throws IOException {
		BlockOutputStream bos = new BlockOutputStream(stream);

		byte[] header = new byte[] { GZIP_MAGIC_LOW, GZIP_MAGIC_HIGH, Deflater.DEFLATED, // CM
				0, // FLG
				0, 0, 0, 0, // MTIME (int)
				0, // XFLG
				-1 // OS
		};

		bos.writeBytes(header);

		return bos;
	}

	@Override
	public void updateChecksumPreCompression(UncompressedBlockData current) {
		if (current == null)
			return;

		byte[] currentData = current.getData();

		if (currentData.length > 0) {
			length += currentData.length;
			crc.update(currentData, 0, currentData.length);
		}
	}

	@Override
	public GzipBlockData compress(UncompressedBlockData previous, UncompressedBlockData current) throws Exception {
		if (current == null)
			return null;

		byte[] rawData = current.getData();

		if (rawData.length == 0)
			return new GzipBlockData(new byte[0]);

		try {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			Deflater def = new Deflater(compressLevel, true);

			if (previous != null)
				def.setDictionary(previous.getData());

			BlockCompressionHelper helper = new BlockCompressionHelper(baos, def);
			helper.write(rawData);

			byte out[] = baos.toByteArray();
			def.finish();

			return new GzipBlockData(out);
		} catch (IOException e) {
			throw new Exception("Gzip compression exception", e);
		}
	}

	@Override
	public void updateChecksumPostCompression(BlockData current) {
	}

	public void writeTrailer(BlockOutputStream stream) throws IOException {
		int crc = (int) this.crc.getValue();
		int length = this.length;

		byte[] trailer = new byte[10];

		trailer[0] = (byte) 0x3; // Empty compressed block with FINAL flag
		trailer[1] = (byte) 0x0;

		trailer[2] = (byte) (crc & 0xFF);
		trailer[3] = (byte) ((crc >> 8) & 0xFF);
		trailer[4] = (byte) ((crc >> 16) & 0xFF);
		trailer[5] = (byte) ((crc >> 24) & 0xFF);

		trailer[6] = (byte) (length & 0xFF);
		trailer[7] = (byte) ((length >> 8) & 0xFF);
		trailer[8] = (byte) ((length >> 16) & 0xFF);
		trailer[9] = (byte) ((length >> 24) & 0xFF);

		stream.writeBytes(trailer);
	}

	private static int BUFFER_SIZE = 512;

	private class BlockCompressionHelper {
		protected OutputStream out;
		protected Deflater def;
		protected byte[] buf;

		public BlockCompressionHelper(OutputStream out, Deflater def) {
			this.out = out;
			this.def = def;
			this.buf = new byte[BUFFER_SIZE];
		}

		public void write(byte[] b) throws IOException {
			def.setInput(b, 0, b.length);

			while (!def.needsInput()) {
				deflate(Deflater.NO_FLUSH);
			}

			int len = deflate(Deflater.SYNC_FLUSH);
			while (len == buf.length)
				len = deflate(Deflater.SYNC_FLUSH);

			out.flush();
		}

		protected int deflate(int mode) throws IOException {
			int len = def.deflate(buf, 0, buf.length, mode);
			if (len > 0) {
				out.write(buf, 0, len);
			}

			return len;
		}

	}

}
