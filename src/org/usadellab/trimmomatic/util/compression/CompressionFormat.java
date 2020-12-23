package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipInputStream;

import org.itadaki.bzip2.BZip2InputStream;
import org.itadaki.bzip2.BZip2OutputStream;
import org.usadellab.trimmomatic.util.Logger;

public enum CompressionFormat
{
	GZIP(".gz"),
	BZIP2(".bz2"),
	ZIP(".zip"),
	NONE("");
	
	private final String ext;
	private static Set<CompressionFormat> extSet;
	
	static 
	{
		extSet=new HashSet<CompressionFormat>();
		
		for (CompressionFormat cf : CompressionFormat.values())
			{
		    if (cf.ext.length()>0)
		    	extSet.add(cf);
			}
	
	}
	
	private CompressionFormat(String ext)
	{
		this.ext = ext;
	}
		
	public static CompressionFormat forName(String name)
	{
		String lcName=name.toLowerCase();
	
		for (CompressionFormat cf : extSet)
			{
			if (lcName.endsWith(cf.ext))
					return cf;
			}
		
		return NONE;
	}
	
	

	public static InputStream wrapStreamForParsing(InputStream is, String name) throws IOException
	{
		CompressionFormat cf = CompressionFormat.forName(name);
		
		switch(cf)
			{
			case GZIP:
				return new ConcatGZIPInputStream(is);
			case BZIP2:
				return new BZip2InputStream(is, false);
			case ZIP:
				return new ZipInputStream(is);
			default:
				return is;
			}
	}

	public static OutputStream wrapStreamForSerializing(OutputStream os, String name, Integer compressLevel) throws IOException
	{
		CompressionFormat cf = CompressionFormat.forName(name);
	
		switch(cf)
			{
			case GZIP:
				if (compressLevel == null)
					return new GZIPOutputStream(os);
				else
					return new TunableGZIPOutputStream(os, compressLevel);
			case BZIP2:
				if (compressLevel == null)
					return new BZip2OutputStream(os); 
				else
					return new BZip2OutputStream(os, compressLevel); 
			default:
				return os;
			}
	}
	
	public static ParallelCompressor parallelCompressorForSerializing(Logger logger, String name, Integer compressLevel)
	{
		CompressionFormat cf = CompressionFormat.forName(name);
	
		switch(cf)
			{
			case GZIP:
				return new GzipParallelCompressor(logger, compressLevel);
			case BZIP2:
				return new Bzip2ParallelCompressor(logger, compressLevel);
			default:
				return null;
			}

	}
	
	

}
