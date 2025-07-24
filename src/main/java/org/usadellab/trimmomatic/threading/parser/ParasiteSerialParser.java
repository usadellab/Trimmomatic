package org.usadellab.trimmomatic.threading.parser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class ParasiteSerialParser extends Parser
{
	public ParasiteSerialParser(FastqParser parser)
	{
		super(parser);	
	}
	
	public List<FastqRecord> poll() throws Exception
	{
		if (isComplete())
			return new ArrayList<FastqRecord>();
		
		return parseBlock();
	}
	

}
