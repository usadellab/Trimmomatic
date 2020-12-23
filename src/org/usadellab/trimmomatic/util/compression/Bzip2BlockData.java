package org.usadellab.trimmomatic.util.compression;

import java.util.List;

public class Bzip2BlockData implements BlockData
{
	private byte[] data;
	private long bitCount;
	private List<Integer> blockCRCs;

	public Bzip2BlockData(byte []data, long bitCount, List<Integer> blockCRCs)
	{
		if (data==null)
			throw new NullPointerException("Bzip2BlockData: cannot be null");
	
		this.data=data;
		this.bitCount=bitCount;
		this.blockCRCs=blockCRCs;
	}

	public byte[] getData()
	{
		return data;
	}
	
	public long getBitCount()
	{
		return bitCount;
	}
	
	public List<Integer> getBlockCRCs()
	{
		return blockCRCs;
	}

}
