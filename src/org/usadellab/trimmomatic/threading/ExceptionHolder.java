package org.usadellab.trimmomatic.threading;

public class ExceptionHolder
{
	private Exception exception;

	public ExceptionHolder()
	{		
	}
	
	public synchronized Exception getException()
	{
		return exception;
	}
	
	public synchronized void rethrow() throws Exception
	{
		if(exception!=null)			
			throw exception;
	}
		
	public synchronized void setException(Exception e)
	{
		if(exception==null)			
			this.exception=e;		
	}
	
}
