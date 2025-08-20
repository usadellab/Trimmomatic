package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TailCropTrimmer implements Trimmer {
	private int firstBases;
	private int midBases;
	private int lastBases;

	private int firstMaxLength = Integer.MAX_VALUE / 2;
	private int midMaxLength = Integer.MAX_VALUE / 2;
	private int lastMaxLength = Integer.MAX_VALUE / 2;

	public TailCropTrimmer(String args) {
		String arg[] = args.split(":");

		switch (arg.length) {
		case 1:
			firstBases = midBases = lastBases = Integer.parseInt(arg[0]);
			break;

		case 2:
			firstBases = midBases = lastBases = Integer.parseInt(arg[0]);
			firstMaxLength = midMaxLength = lastMaxLength = Integer.parseInt(arg[1]);
			break;

		case 4:
			firstBases = midBases = Integer.parseInt(arg[0]);
			firstMaxLength = midMaxLength = Integer.parseInt(arg[1]);
			lastBases = Integer.parseInt(arg[2]);
			lastMaxLength = Integer.parseInt(arg[3]);
			break;

		case 6:
			firstBases = Integer.parseInt(arg[0]);
			firstMaxLength = Integer.parseInt(arg[1]);
			midBases = Integer.parseInt(arg[2]);
			midMaxLength = Integer.parseInt(arg[3]);
			lastBases = Integer.parseInt(arg[4]);
			lastMaxLength = Integer.parseInt(arg[5]);
			break;

		default:

		}
	}

	@Override
	public FastqRecord[] processRecords(FastqRecord[] in) {

		if (in == null)
			return null;

		FastqRecord out[] = new FastqRecord[in.length];

		int lastRecord = in.length - 1;

		if (in.length > 0) {
			if (in[0] != null)
				out[0] = processRecord(in[0], firstBases, firstMaxLength);
		}

		for (int i = 1; i < lastRecord; i++) {
			if (in[i] != null)
				out[i] = processRecord(in[i], midBases, midMaxLength);
		}

		if (lastRecord > 0) {
			if (in[lastRecord] != null)
				out[lastRecord] = processRecord(in[lastRecord], lastBases, lastMaxLength);
		}

		return out;
	}

	private FastqRecord processRecord(FastqRecord in, int bases, int maxLength) {
		int len = in.getSequence().length();

		int toTrim = bases;
		int overLen = len - toTrim - maxLength;

		if (overLen > 0)
			toTrim += overLen;

		if (len <= toTrim)
			return null;

		if (toTrim == 0)
			return in;

		return new FastqRecord(in, 0, len - toTrim);
	}

}
