package org.usadellab.trimmomatic.fastq.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public interface Trimmer
{
	public FastqRecord[] processRecords(FastqRecord in[]);
}
