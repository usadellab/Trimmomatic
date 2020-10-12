package org.usadellab.trimmomatic.util;

public class Logger
{
	private boolean showError, showWarning, showInfo;

	public Logger(boolean showError, boolean showWarning, boolean showInfo)
	{
		this.showError=showError;
		this.showWarning=showWarning;
		this.showInfo=showInfo;
	}

	public void handleException(Exception ex)
	{
		ex.printStackTrace(System.err);
	}
	
	public void errorln()
	{
		if(showError)		
			System.err.println();
	}
	
	public void errorln(String message)
	{
		if(showError)		
			System.err.println(message);
	}

	public void error(String message)
	{
		if(showError)		
			System.err.print(message);
	}
	
	

	
	
	public void warnln()
	{
		if(showWarning)		
			System.err.println();
	}
	
	public void warnln(String message)
	{
		if(showWarning)		
			System.err.println(message);
	}

	public void warn(String message)
	{
		if(showWarning)		
			System.err.print(message);
	}

	
	public void infoln()
	{
		if(showInfo)
			System.err.println();
	}
	
	public void infoln(String message)
	{
		if(showInfo)
			System.err.println(message);
	}

	public void info(String message)
	{
		if(showInfo)
			System.err.print(message);
	}

}
