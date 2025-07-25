package org.usadellab.trimmomatic.util.compression;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.itadaki.bzip2.BZip2BlockCompressor;
import org.itadaki.bzip2.BitOutputStream;
import org.usadellab.trimmomatic.util.Logger;

public class Bzip2ParallelCompressor implements ParallelCompressor
{
	private static final int DEFAULT_COMPRESSION_LEVEL = 9;
	private static final int BLOCK_SIZE_MULTIPLER = 100000;

	private static final int HEADER_MAGIC_1=0x425a;
	private static final int HEADER_MAGIC_2=0x68;
	
	private static final int TRAILER_MAGIC_1 = 0x177245;
	private static final int TRAILER_MAGIC_2 = 0x385090;
	
	@SuppressWarnings("unused")
	private Logger logger;	
	private int compressLevel;
	private int block_size;
	
	private int streamCRC;
	
	public Bzip2ParallelCompressor(Logger logger, Integer compressLevel)
	{
		this.logger = logger;
		this.compressLevel = compressLevel==null ? DEFAULT_COMPRESSION_LEVEL : compressLevel;
		this.block_size = this.compressLevel * BLOCK_SIZE_MULTIPLER;
				
		streamCRC=0;
	}
	
	
	
	@Override
	public BlockOutputStream wrapAndWriteHeader(OutputStream stream) throws IOException
	{
		Bzip2BlockOutputStream bos = new Bzip2BlockOutputStream(stream);
		
		bos.writeBits (16, HEADER_MAGIC_1);
		bos.writeBits (8,  HEADER_MAGIC_2);
		bos.writeBits (8, '0' + compressLevel);
		
		return bos;
	}

	@Override
	public void updateChecksumPreCompression(UncompressedBlockData data)
	{
	
	}

	@Override
	public Bzip2BlockData compress(UncompressedBlockData previous, UncompressedBlockData current) throws Exception
	{
		if (current==null)
			return null;
	
		byte[] rawData = current.getData();
	
		if(rawData.length==0)	
			return new Bzip2BlockData(new byte[0], 0, new ArrayList<Integer>());
	
		try
			{		
			BlockBitOutputStream baos=new BlockBitOutputStream();
			BlockCompressionHelper helper = new BlockCompressionHelper(baos);
		
			helper.write(rawData);					
			List<Integer> blockCRCs = helper.getBlockCRC();
			
			long bitCount = baos.getBitCount();
			baos.flush();
			
			byte out[]=baos.getOutputStream().toByteArray();			
			
			return new Bzip2BlockData(out, bitCount, blockCRCs);
			}		
		catch (IOException e)
			{
			throw new Exception("Bzip2 compression exception", e);
			}
	}
	
	@Override
	public void updateChecksumPostCompression(BlockData data)
	{
		Bzip2BlockData blockData = (Bzip2BlockData) data;
		
		List<Integer> blockCRCs=blockData.getBlockCRCs();
		
		for(int blockCRC: blockCRCs)		
			streamCRC = ((this.streamCRC << 1) | (this.streamCRC >>> 31)) ^ blockCRC;		
	}
	
	@Override
	public void writeTrailer(BlockOutputStream stream) throws IOException
	{
		Bzip2BlockOutputStream bbos = (Bzip2BlockOutputStream) stream;
		
		bbos.writeBits(24, TRAILER_MAGIC_1);
		bbos.writeBits(24, TRAILER_MAGIC_2);
		bbos.writeInteger(this.streamCRC);	
	}

	
	
	private static class Bzip2BlockOutputStream extends BlockOutputStream
	{
		private BufferedOutputStream bufferedStream;
		private BitOutputStream bitOutputStream;
	
		Bzip2BlockOutputStream(OutputStream stream)
		{
			super(stream);
			
			bufferedStream = new BufferedOutputStream(stream, 4096);
			bitOutputStream=new BitOutputStream(bufferedStream);
		}

		private void writeBits(int count, int value) throws IOException
		{
			bitOutputStream.writeBits(count, value);
		}
		
		public void writeInteger (final int value) throws IOException 
		{
			writeBits (16, (value >>> 16) & 0xffff);
			writeBits (16, value & 0xffff);
		}		
		
		@Override
		public void writeBlock(BlockData blockData) throws IOException
		{
			Bzip2BlockData bBlockData = (Bzip2BlockData) blockData;		
			byte[] data = bBlockData.getData();
			long bitCount = bBlockData.getBitCount();
			
			int offset = 0;
			
			while (bitCount>7)
				{
				bitOutputStream.writeBits(8, data[offset]);
				bitCount-=8;
				offset++;
				}
			
			if (bitCount > 0)
				bitOutputStream.writeBits((int)bitCount, data[offset] >>> (8-bitCount) ); // Bits MSB->LSB
		}

		@Override
		public void close() throws IOException
		{
			bitOutputStream.flush();
			bufferedStream.close();			
		}		
		
	}

	
	private static class BlockBitOutputStream extends BitOutputStream
	{
		private ByteArrayOutputStream os;
		private long bitCount;
		
		BlockBitOutputStream()
		{
			this(new ByteArrayOutputStream());
		}
		
		BlockBitOutputStream(ByteArrayOutputStream stream)
		{
			super(stream);
			os=stream;
			bitCount=0;
		}
		
		private ByteArrayOutputStream getOutputStream()
		{
			return os;
		}
		
		private long getBitCount()
		{
			return bitCount;
		}
		
		
		@Override
		public void flush() throws IOException
		{
			super.flush();
		}

		@Override
		public void writeBits(int count, int value) throws IOException
		{
			bitCount+=count;
			super.writeBits(count, value);
		}

		@Override
		public void writeBoolean(boolean value) throws IOException
		{
			bitCount+=1;
			super.writeBoolean(value);
		}
	
	}
	
	
	private class BlockCompressionHelper
	{
		private BlockBitOutputStream bitOutputStream;
		private BZip2BlockCompressor blockCompressor;
		
		private ArrayList<Integer> blockCRCs;
		
		BlockCompressionHelper(BlockBitOutputStream stream) throws IOException
		{
			bitOutputStream = stream;
			blockCRCs=new ArrayList<Integer>();	
		}
		
		public List<Integer> getBlockCRC()
		{
			return blockCRCs;
		}
		
        public void write(byte[] b) throws IOException 
        {
        	initBlock();
        
        	int offset = 0;
        	int length = b.length;
        	
        	int bytesWritten;
        	while (length > 0) 
        		{
				if ((bytesWritten = this.blockCompressor.write (b, offset, length)) < length) 
					{
					closeBlock();
					initBlock();
					}
				offset += bytesWritten;
				length -= bytesWritten;
        		}
        	
        	closeBlock();
        }
		
        private void initBlock() throws IOException
        {
        	blockCompressor = new BZip2BlockCompressor(bitOutputStream, block_size);	
        }
        
    	private void closeBlock() throws IOException 
    	{
    		if (this.blockCompressor.isEmpty()) 
    			return;    		

    		this.blockCompressor.close();    		
    		blockCRCs.add(this.blockCompressor.getCRC());
    	}
        
	}

}
