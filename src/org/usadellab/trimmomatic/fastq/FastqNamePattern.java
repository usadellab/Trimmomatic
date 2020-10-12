package org.usadellab.trimmomatic.fastq;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/*

@HWUSI-EAS100R:6:73:941:1973#0/1
  
HWUSI-EAS100R	the unique instrument name
6	flowcell lane
73	tile number within the flowcell lane
941	'x'-coordinate of the cluster within the tile
1973	'y'-coordinate of the cluster within the tile
#0	index number for a multiplexed sample (0 for no indexing)
/1	the member of a pair, /1 or /2 (paired-end or mate-pair reads only)


EAS139	the unique instrument name
136	the run id
FC706VJ	the flowcell id
2	flowcell lane
2104	tile number within the flowcell lane
15343	'x'-coordinate of the cluster within the tile
197393	'y'-coordinate of the cluster within the tile
1	the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
Y	Y if the read is filtered, N otherwise
18	0 when none of the control bits are on, otherwise it is an even number
ATCACG	index sequence

Also allow for SRA padding:

@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36

*/

// CASAVA_13("[@ ]([^:]*:[0-9]*:[0-9]*:[0-9]*:[0-9]*#[0-9]*)/[12]","$1"), // HWUSI-EAS100R:6:73:941:1973#0/1 
// SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36

public enum FastqNamePattern
{
	CASAVA_13("(.* )?([^:]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+#[A-Z0-9]+).*","$2"), // HWUSI-EAS100R:6:73:941:1973#0/1 -> HWUSI-EAS100R:6:73:941:1973 or HWUSI-EAS100R:6:73:941:1973#NNNN/1 
	CASAVA_18("(.* )?([^:]+:[0-9]+:[A-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+).*","$2");  // EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	
	private Pattern pattern;
	private String replacement;
	
	FastqNamePattern(String patternStr, String replacement)
	{
		this.pattern=Pattern.compile(patternStr);
		this.replacement=replacement;
	}
	
	public boolean match(String str)
	{
		return pattern.matcher(str).find();
	}
	
	public String canonicalizeOne(String str)
	{
		Matcher matcher=pattern.matcher(str);
		if(!matcher.find())
			return null;
		
		return matcher.replaceAll(replacement);
			
	}
	
	public static String canonicalize(String str)
	{
		for(FastqNamePattern fnp: FastqNamePattern.values())
		{
			String canon=fnp.canonicalizeOne(str);

			if(canon!=null)
				return canon;
		}
	
		return null;
	}
	
}
