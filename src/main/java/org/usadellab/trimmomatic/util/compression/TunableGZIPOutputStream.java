package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.GZIPOutputStream;

public class TunableGZIPOutputStream extends GZIPOutputStream {
	public TunableGZIPOutputStream(OutputStream out, int compressLevel) throws IOException {
		super(out);
		def.setLevel(compressLevel);
	}
}
