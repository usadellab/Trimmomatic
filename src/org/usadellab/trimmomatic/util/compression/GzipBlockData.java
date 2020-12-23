package org.usadellab.trimmomatic.util.compression;

public class GzipBlockData implements BlockData
{
	private byte[] data;

	public GzipBlockData(byte[] data)
	{
		if (data==null)
			throw new NullPointerException("GzipBlockData: cannot be null");
	
		this.data=data;
	}

	public byte[] getData()
	{
		return data;
	}

}
