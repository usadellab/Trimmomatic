package org.usadellab.trimmomatic.fasta;

public class FastaRecord
{
	private String name;
	private String sequence;
	private String fullName;
        private String barcodeLabel = null;
	
	public FastaRecord(String name, String sequence)
	{
		this.name=name;
		this.sequence=sequence;                
	}
	
	public FastaRecord(String name, String sequence, String fullName)
	{
		this (name, sequence );
		this.fullName = fullName;
	}

	public String getName()
	{
		return name;
	}

	public String getSequence()
	{
		return sequence;
	}
	
	public String getFullName(){
		return fullName;
	}

        public String getBarcodeLabel() {
            return barcodeLabel;
        }

        public void setBarcodeLabel(String barcodeLabel) {
            this.barcodeLabel = barcodeLabel;
        }
        
	public String getComplementSequence()
	{
		StringBuilder sb=new StringBuilder(sequence);
		
		int len=sequence.length();
		int last=len-1;
		
		for(int i=0;i<sequence.length();i++)
			sb.setCharAt(i, complementChar(sequence.charAt(last-i)));
		
		return sb.toString();
	}
	
	
	
	public String getEndSubsequence(int length, boolean endOfRead, boolean complement)
	{
		String sub;
		if(endOfRead)
			sub=sequence.substring(sequence.length()-length);
		else
			sub=sequence.substring(0,length);
		
		if(complement)
			{
			StringBuilder sb=new StringBuilder(sub);
	
			int len=sub.length();
			int last=len-1;
			
			for(int i=0;i<len;i++)
				sb.setCharAt(i, complementChar(sub.charAt(last-i)));
			
			sub=sb.toString();
			}
	
		return sub;
	}

	
	public static char complementChar(char ch)
	{
		switch (ch)
			{
			case 'A':
				return 'T';
			case 'C':
				return 'G';
			case 'G':
				return 'C';
			case 'T':
				return 'A';
			}

		return 'N';
	}
	
	
	public static String getComplementSequence(String sequence)
	{
		StringBuilder sb=new StringBuilder(sequence);
		
		int len=sequence.length();
		int last=len-1;
		
		for(int i=0;i<sequence.length();i++)
			sb.setCharAt(i, complementChar(sequence.charAt(last-i)));
		
		return sb.toString();
	}
	
}
