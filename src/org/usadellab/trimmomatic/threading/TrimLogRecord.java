/**
 * 
 */
package org.usadellab.trimmomatic.threading;

public class TrimLogRecord
{
	private String readName;
	private int length;
	private int startPos;
	private int endPos;
	private int trimTail;
	
	public TrimLogRecord(String readName, int length, int startPos, int endPos, int trimTail)
	{
		this.readName = readName;
		this.length = length;
		this.startPos = startPos;
		this.endPos = endPos;
		this.trimTail = trimTail;
	}
	public String getReadName()
	{
		return readName;
	}
	public int getLength()
	{
		return length;
	}
	public int getStartPos()
	{
		return startPos;
	}
	public int getEndPos()
	{
		return endPos;
	}
	public int getTrimTail()
	{
		return trimTail;
	}
}