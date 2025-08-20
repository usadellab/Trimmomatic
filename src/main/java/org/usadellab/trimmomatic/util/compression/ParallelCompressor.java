package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.OutputStream;

public interface ParallelCompressor {
	BlockOutputStream wrapAndWriteHeader(OutputStream stream) throws IOException;

	void updateChecksumPreCompression(UncompressedBlockData current);

	BlockData compress(UncompressedBlockData previous, UncompressedBlockData current) throws Exception;

	void updateChecksumPostCompression(BlockData current);

	void writeTrailer(BlockOutputStream stream) throws IOException;
}
