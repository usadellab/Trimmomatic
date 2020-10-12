package org.usadellab.trimmomatic.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.PushbackInputStream;
import java.util.zip.GZIPInputStream;

public class ConcatGZIPInputStream extends InputStream
{
	private PushbackInputStream source;
	private GZIPHelperInputStream gzIn;

	public ConcatGZIPInputStream(InputStream in) throws IOException
	{
		source = new PushbackInputStream(in, 1024);
		nextGzipInputStream();
	}

	private void nextGzipInputStream() throws IOException
	{
		boolean more=false;
		
		if((gzIn!=null)&&(gzIn.pushbackUnused()>0))
				more=true;
	
		if(!more)
			{
			int r=source.read();
			if(r!=-1)
				{
				source.unread(r);
				more=true;
				}
			}
		
		if(more)
			{
			gzIn=new GZIPHelperInputStream(source);
			}
		else
			gzIn=null;
	}
		
	@Override
	public void close() throws IOException
	{
		gzIn=null;
		source.close();
	}

	@Override
	public int read() throws IOException
	{
		int res=-1;
		
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read();
			if(res==-1)
				nextGzipInputStream();
			}

		return res;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException
	{
		int res=-1;
	
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read(b,off,len);
			if(res==-1)
				nextGzipInputStream();
			}
		
		return res;
	}

	@Override
	public int read(byte[] b) throws IOException
	{
		int res=-1;
	
		while(res==-1 && gzIn!=null)
			{
			res=gzIn.read(b);
			if(res==-1)
				nextGzipInputStream();
			}

		return res;	
	}

	private class GZIPHelperInputStream extends GZIPInputStream
	{
		private GZIPHelperInputStream(InputStream in) throws IOException
		{
			super(in);
		}
		
		private int pushbackUnused() throws IOException
		{
			int amount=inf.getRemaining()-8;
			if(amount>0)
				source.unread(buf, len-amount, amount);
			
			return amount;
		}
	}

}
