package org.usadellab.trimmomatic.trim;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.usadellab.trimmomatic.fasta.FastaParser;
import org.usadellab.trimmomatic.fasta.FastaRecord;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmer implements Trimmer
{
	public final static String PREFIX = "Prefix";
	public final static String SUFFIX_F = "/1";
	public final static String SUFFIX_R = "/2";

	public final static int INTERLEAVE = 4;

	private final static float LOG10_4 = 0.60206f;

	private Logger logger;
	
	private int seedMaxMiss;
	private int minPalindromeLikelihood;
	// private int minPalindromeOverlap;

	private int minSequenceLikelihood;
	private int minSequenceOverlap;

	private int minPrefix;
	public boolean palindromeKeepBoth;

	private List<IlluminaPrefixPair> prefixPairs;

	private Set<IlluminaClippingSeq> forwardSeqs;
	private Set<IlluminaClippingSeq> reverseSeqs;
	private Set<IlluminaClippingSeq> commonSeqs;

	
	public static IlluminaClippingTrimmer makeIlluminaClippingTrimmer(Logger logger, String args) throws IOException
	{
		String arg[] = args.split(":");
		
		File seqs=new File(arg[0]);

		int seedMaxMiss = Integer.parseInt(arg[1]);
		int minPalindromeLikelihood = Integer.parseInt(arg[2]);
		int minSequenceLikelihood = Integer.parseInt(arg[3]);

		int minPrefix=8;
		boolean palindromeKeepBoth=false;
		
		if (arg.length > 4)
			minPrefix = Integer.parseInt(arg[4]);

		if (arg.length > 5)
			palindromeKeepBoth = Boolean.parseBoolean(arg[5]);

		IlluminaClippingTrimmer trimmer=new IlluminaClippingTrimmer(logger, seedMaxMiss, minPalindromeLikelihood, minSequenceLikelihood, minPrefix, palindromeKeepBoth);

		try
			{
			trimmer.loadSequences(seqs.getCanonicalPath());
			}
		catch (IOException ex)
			{
			logger.handleException(ex);
			}
		
		return trimmer;
	}


	public IlluminaClippingTrimmer(Logger logger, int seedMaxMiss, int minPalindromeLikelihood, int minSequenceLikelihood, int minPrefix, boolean palindromeKeepBoth)
	{
		this.logger=logger;
		
		this.seedMaxMiss = seedMaxMiss;
		this.minPalindromeLikelihood = minPalindromeLikelihood;
		this.minSequenceLikelihood = minSequenceLikelihood;
		
		this.minPrefix=minPrefix;
		this.palindromeKeepBoth=palindromeKeepBoth;
		
		// minPalindromeOverlap=(int)(minPalindromeLikelihood/LOG10_4);
		minSequenceOverlap = (int) (minSequenceLikelihood / LOG10_4);

		if (minSequenceOverlap > 15)
			minSequenceOverlap = 15;
		
		prefixPairs = new ArrayList<IlluminaPrefixPair>();
		commonSeqs = new HashSet<IlluminaClippingSeq>();
		forwardSeqs = new HashSet<IlluminaClippingSeq>();
		reverseSeqs = new HashSet<IlluminaClippingSeq>();		
	}

	
	
	private void loadSequences(String seqPath) throws IOException
	{
		FastaParser parser = new FastaParser();
		parser.parse(new File(seqPath));

		Map<String, FastaRecord> forwardSeqMap = new HashMap<String, FastaRecord>();
		Map<String, FastaRecord> reverseSeqMap = new HashMap<String, FastaRecord>();
		Map<String, FastaRecord> commonSeqMap = new HashMap<String, FastaRecord>();

		Set<String> forwardPrefix = new HashSet<String>();
		Set<String> reversePrefix = new HashSet<String>();

		while (parser.hasNext())
			{
			FastaRecord rec = parser.next();

			String name = rec.getName();

			if (name.endsWith(SUFFIX_F))
				{
				forwardSeqMap.put(name, rec);

				if (name.startsWith(PREFIX))
					{
					String clippedName = name.substring(0, name.length() - SUFFIX_F.length());
					forwardPrefix.add(clippedName);
					}
				}
			else if (name.endsWith(SUFFIX_R))
				{
				reverseSeqMap.put(name, rec);

				if (name.startsWith(PREFIX))
					{
					String clippedName = name.substring(0, name.length() - SUFFIX_R.length());
					reversePrefix.add(clippedName);
					}
				}
			else
				commonSeqMap.put(name, rec);
			}

		Set<String> prefixSet = new HashSet<String>(forwardPrefix);
		prefixSet.retainAll(reversePrefix);

		prefixPairs = new ArrayList<IlluminaPrefixPair>();

		for (String prefix : prefixSet)
			{
			String forwardName = prefix + SUFFIX_F;
			String reverseName = prefix + SUFFIX_R;

			FastaRecord forwardRec = forwardSeqMap.remove(forwardName);
			FastaRecord reverseRec = reverseSeqMap.remove(reverseName);

			prefixPairs.add(new IlluminaPrefixPair(forwardRec.getSequence(), reverseRec.getSequence()));
			}

		forwardSeqs = mapClippingSet(forwardSeqMap);
		reverseSeqs = mapClippingSet(reverseSeqMap);
		commonSeqs = mapClippingSet(commonSeqMap);

		logger.infoln("ILLUMINACLIP: Using " + prefixPairs.size() + " prefix pairs, " + commonSeqs.size()
				+ " forward/reverse sequences, " + forwardSeqs.size() + " forward only sequences, "
				+ reverseSeqs.size() + " reverse only sequences");

	}

	
	void addPrefixPair(String prefix1, String prefix2)
	{
		prefixPairs.add(new IlluminaPrefixPair(prefix1, prefix2));
	}

	void addClippingSeq(IlluminaClippingSeq clippingSeq, boolean forward, boolean reverse)
		{
		if(forward)
			{
			if(reverse)
				commonSeqs.add(clippingSeq);
			else
				forwardSeqs.add(clippingSeq);
			}
		else
			reverseSeqs.add(clippingSeq);
		}
	
	
	
	private Set<IlluminaClippingSeq> mapClippingSet(Map<String, FastaRecord> map)
	{
		Set<String> uniqueSeq = new HashSet<String>();
		Set<IlluminaClippingSeq> out = new HashSet<IlluminaClippingSeq>();

		for (FastaRecord rec : map.values())
			{
			String seq = rec.getSequence();

			if (uniqueSeq.contains(seq))
				logger.warnln("Skipping duplicate Clipping Sequence: '" + seq + "'");
			else
				{
				uniqueSeq.add(seq);

				if (seq.length() < 16)
					out.add(new IlluminaShortClippingSeq(rec.getSequence()));
				else if (seq.length() < 24)
					out.add(new IlluminaMediumClippingSeq(rec.getSequence()));
				else
					out.add(new IlluminaLongClippingSeq(rec.getSequence()));
				}

			}

		return out;
	}

	private Integer min(Integer a, Integer b)
	{
		if (a == null)
			return b;

		if (b == null)
			return a;

		return (a < b) ? a : b;
	}

	@Override
	public FastqRecord[] processRecords(FastqRecord[] in)
	{
		FastqRecord forwardRec = null;
		FastqRecord reverseRec = null;

		if (in.length > 0)
			forwardRec = in[0];

		if (in.length > 1)
			reverseRec = in[1];

		Integer toKeepForward = null, toKeepReverse = null;

		if (forwardRec != null && reverseRec != null)
			{
			// First, check for a palindrome

			for (IlluminaPrefixPair pair : prefixPairs)
				{
				Integer toKeep = pair.palindromeReadsCompare(forwardRec, reverseRec);

				if (toKeep != null)
					{
					toKeepForward = min(toKeepForward, toKeep);

					if (palindromeKeepBoth)
						toKeepReverse = min(toKeepReverse, toKeep);
					else
						toKeepReverse = 0;
					}
				}
			}

		// Also check each record for other nasties

		if (forwardRec != null)
			{
			if (toKeepForward == null || toKeepForward > 0)
				{
				for (IlluminaClippingSeq seq : forwardSeqs)
					toKeepForward = min(toKeepForward, seq.readsSeqCompare(forwardRec));

				for (IlluminaClippingSeq seq : commonSeqs)
					toKeepForward = min(toKeepForward, seq.readsSeqCompare(forwardRec));
				}

			// Keep the minimum

			if (toKeepForward != null)
				{
				if (toKeepForward > 0)
					forwardRec = new FastqRecord(forwardRec, 0, toKeepForward);
				else
					forwardRec = null;
				}
			}

		if (reverseRec != null)
			{
			if (toKeepReverse == null || toKeepReverse > 0)
				{
				for (IlluminaClippingSeq seq : reverseSeqs)
					toKeepReverse = min(toKeepReverse, seq.readsSeqCompare(reverseRec));

				for (IlluminaClippingSeq seq : commonSeqs)
					toKeepReverse = min(toKeepReverse, seq.readsSeqCompare(reverseRec));
				}
			// Keep the minimum

			if (toKeepReverse != null)
				{
				if (toKeepReverse > 0)
					reverseRec = new FastqRecord(reverseRec, 0, toKeepReverse);
				else
					reverseRec = null;
				}
			}

		if (in.length == 2)
			return new FastqRecord[] { forwardRec, reverseRec };
		if (in.length == 1)
			return new FastqRecord[] { forwardRec };

		return new FastqRecord[0];
	}

	class IlluminaPrefixPair
	{
		private String prefix1;
		private String prefix2;

		private IlluminaPrefixPair(String prefix1, String prefix2)
		{
			logger.infoln("Using PrefixPair: '" + prefix1 + "' and '" + prefix2 + "'");

			int length1 = prefix1.length();
			int length2 = prefix2.length();

			if (length1 != length2)
				{
				int minLength = length1;
				if (length2 < minLength)
					minLength = length2;

				prefix1 = prefix1.substring(length1 - minLength);
				prefix2 = prefix2.substring(length2 - minLength);
				}

			this.prefix1 = prefix1;
			this.prefix2 = prefix2;
		}

		public String getPrefix1()
		{
			return prefix1;
		}

		public String getPrefix2()
		{
			return prefix2;
		}


		
		private Integer palindromeReadsCompare(FastqRecord rec1, FastqRecord rec2)
		{
			int seedMax = seedMaxMiss * 2;

			long pack1[] = packSeqInternal(getPrefix1() + rec1.getSequence(), false);
			long pack2[] = packSeqInternal(getPrefix2() + rec2.getSequence(), true);

			int prefixLength = getPrefix1().length();

			int testIndex = 0;
			int refIndex = prefixLength;

			if(pack1.length<=refIndex || pack2.length<=refIndex)
				return null;
			
			int count = 0;

			int seedSkip = prefixLength - 16;
			if (seedSkip > 0)
				{
				testIndex = seedSkip;
				count = seedSkip;
				}

			long ref1, ref2;

			int seqlen1 = rec1.getSequence().length() + prefixLength;
			int seqlen2 = rec2.getSequence().length() + prefixLength;

			int maxCount = (seqlen1 > seqlen2 ? seqlen1 : seqlen2) - 15 - minPrefix;

			while (count < maxCount)
				{
				ref1 = pack1[refIndex];
				ref2 = pack2[refIndex];

				if ((testIndex < pack2.length && Long.bitCount(ref1 ^ pack2[testIndex]) <= seedMax)
						|| (testIndex < pack1.length && Long.bitCount(ref2 ^ pack1[testIndex]) <= seedMax))				
					{
					int totalOverlap = count + prefixLength + 16;

					int skip1 = 0;
					int skip2 = 0;

					if (totalOverlap > seqlen1)
						skip2 = totalOverlap - seqlen1;

					if (totalOverlap > seqlen2)
						skip1 = totalOverlap - seqlen2;

					int actualOverlap = totalOverlap - skip1 - skip2;

					float palindromeLikelihood = calculatePalindromeDifferenceQuality(rec1, rec2, actualOverlap, skip1,
							skip2);

					if (palindromeLikelihood >= minPalindromeLikelihood)
						{
						return totalOverlap - prefixLength * 2;
						}
					}

				count++;
				int testRefIndex = refIndex + 1;

				if (((count & 0x1) == 0) && testRefIndex < pack1.length && testRefIndex < pack2.length)
					refIndex++;
				else
					testIndex++;

				}

			return null;
		}

		private char compCh(char ch)
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

		private float calculatePalindromeDifferenceQuality(FastqRecord rec1, FastqRecord rec2, int overlap, int skip1,
				int skip2)
		{
			String seq1 = rec1.getSequence();
			String seq2 = rec2.getSequence();

			String prefix1 = getPrefix1();
			String prefix2 = getPrefix2();

			int quals1[] = rec1.getQualityAsInteger(true);
			int quals2[] = rec2.getQualityAsInteger(true);

			int prefixLength = prefix1.length();

			float likelihood[] = new float[overlap];

			for (int i = 0; i < overlap; i++)
				{
				int offset1 = i + skip1;
				int offset2 = skip2 + overlap - i - 1;

				char ch1 = offset1 < prefixLength ? prefix1.charAt(offset1) : seq1.charAt(offset1 - prefixLength);
				char ch2 = offset2 < prefixLength ? prefix2.charAt(offset2) : seq2.charAt(offset2 - prefixLength);

				ch2 = compCh(ch2);

				int qual1 = offset1 < prefixLength ? 100 : quals1[offset1 - prefixLength];
				int qual2 = offset2 < prefixLength ? 100 : quals2[offset2 - prefixLength];

				if (ch1 == 'N' || ch2 == 'N')
					likelihood[i] = 0;
				else if (ch1 != ch2)
					{
					if (qual1 < qual2)
						likelihood[i] = -qual1 / 10;
					else
						likelihood[i] = -qual2 / 10;
					}
				else
					likelihood[i] = LOG10_4;
				}

			return calculateTotal(likelihood);
		}

		private float calculateTotal(float vals[])
		{
			float total = 0;

			for (float val : vals)
				total += val;

			return total;
		}
	}

	abstract class IlluminaClippingSeq
	{
		String seq;
		long pack[];

		
		public String getSeq()
		{
			return seq;
		}

		public long[] getPack()
		{
			return pack;
		}

		abstract Integer readsSeqCompare(FastqRecord rec);
		

		float calculateDifferenceQuality(FastqRecord rec, String clipSeq, int overlap, int recOffset)
		{
			String seq = rec.getSequence();
			int quals[] = rec.getQualityAsInteger(true);

			int recPos = (recOffset > 0) ? recOffset : 0;
			int clipPos = (recOffset < 0) ? -recOffset : 0;

			float likelihood[] = new float[overlap];

			for (int i = 0; i < overlap; i++)
				{
				char ch1 = seq.charAt(recPos);
				char ch2 = clipSeq.charAt(clipPos);

				if ((ch1 == 'N') || (ch2 == 'N'))
					{
					likelihood[i] = 0;
					}
				else if (ch1 != ch2)
					{
					likelihood[i] = -quals[recPos] / 10.0f;
					}
				else
					likelihood[i] = LOG10_4;

				recPos++;
				clipPos++;
				}

			float l = calculateMaximumRange(likelihood);

			return l;
		}

		private float calculateMaximumRange(float vals[])
		{
			List<Float> merges = new ArrayList<Float>();
			float total = 0;

			for (float val : vals)
				{
				if ((total > 0 && val < 0) || (total < 0 && val > 0))
					{
					merges.add(total);
					total = val;
					}
				else
					total += val;
				}
			merges.add(total);

			boolean scanAgain = true;

			while (merges.size() > 0 && scanAgain)
				{
				ListIterator<Float> mergeIter = merges.listIterator();
				scanAgain = false;

				while (mergeIter.hasNext())
					{
					float val = mergeIter.next();

					if (val < 0 && mergeIter.hasPrevious() && mergeIter.hasNext())
						{
						float prev = mergeIter.previous();
						mergeIter.next();
						float next = mergeIter.next();

						if ((prev > -val) && (next > -val))
							{
							mergeIter.remove();
							mergeIter.previous();
							mergeIter.remove();
							mergeIter.previous();
							mergeIter.set(prev + val + next);

							scanAgain = true;
							}
						else
							mergeIter.previous();
						}
					}

				}

			float max = 0;
			for (float val : merges)
				{
				if (val > max)
					max = val;
				}

			return max;
		}
	}

	class IlluminaShortClippingSeq extends IlluminaClippingSeq
	{
		private long mask;

		IlluminaShortClippingSeq(String seq)
		{
			logger.infoln("Using Short Clipping Sequence: '" + seq + "'");

			this.seq = seq;
			this.mask = calcSingleMask(seq.length());
			pack = packSeqExternal(seq);
		}

		public long getMask()
		{
			return mask;
		}

		public Integer readsSeqCompare(FastqRecord rec)
		{
			int seedMax = seedMaxMiss * 2;

			String recSequence = rec.getSequence();
			String clipSequence = seq;

			Set<Integer> offsetSet = new TreeSet<Integer>();

			long packRec[] = packSeqExternal(rec.getSequence());
			long packClip[] = getPack();
			long mask = getMask();

			int packRecMax = packRec.length - minSequenceOverlap;
			int packClipMax = packClip.length - minSequenceOverlap;

			for (int i = 0; i < packRecMax; i++)
				{
				long comboMask = calcSingleMask(packRec.length - i) & mask;

				for (int j = 0; j < packClipMax; j++)
					{
					int diff = Long.bitCount((packRec[i] ^ packClip[j]) & comboMask);

					if (diff <= seedMax)
						{
						int offset = i - j;
						offsetSet.add(offset);
						}
					}
				}

			for (Integer offset : offsetSet)
				{
				int recCompLength = offset > 0 ? recSequence.length() - offset : recSequence.length();
				int clipCompLength = offset < 0 ? clipSequence.length() + offset : clipSequence.length();

				int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

				if (compLength > minSequenceOverlap)
					{
					float seqLikelihood = calculateDifferenceQuality(rec, clipSequence, compLength, offset);

					if (seqLikelihood >= minSequenceLikelihood)
						return offset;
					}
				}

			return null;
		}
	}

	class IlluminaMediumClippingSeq extends IlluminaClippingSeq
	{
		IlluminaMediumClippingSeq(String seq)
		{
			logger.infoln("Using Medium Clipping Sequence: '" + seq + "'");

			this.seq = seq;
			pack = packSeqInternal(seq, false);
		}
		
		
		public Integer readsSeqCompare(FastqRecord rec)
		{
			int seedMax = seedMaxMiss * 2;

			String recSequence = rec.getSequence();
			String clipSequence = seq;

			Set<Integer> offsetSet = new TreeSet<Integer>();

			long packRec[] = packSeqExternal(rec.getSequence());
			long packClip[] = getPack();

			int packRecMax = packRec.length - minSequenceOverlap;
			int packClipMax = packClip.length;

			for (int i = 0; i < packRecMax; i++)
				{
				long comboMask = calcSingleMask(packRec.length - i);

				for (int j = 0; j < packClipMax; j++)
					{
					int diff = Long.bitCount((packRec[i] ^ packClip[j]) & comboMask);

					if (diff <= seedMax)
						{
						int offset = i - j;
						offsetSet.add(offset);
						}
					}
				}

			for (Integer offset : offsetSet)
				{
				int recCompLength = offset > 0 ? recSequence.length() - offset : recSequence.length();
				int clipCompLength = offset < 0 ? clipSequence.length() + offset : clipSequence.length();

				int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

				if (compLength > minSequenceOverlap)
					{
					float seqLikelihood = calculateDifferenceQuality(rec, clipSequence, compLength, offset);

					if (seqLikelihood >= minSequenceLikelihood)
						return offset;
					}
				}

			return null;
		}
		

	}
	
	
	class IlluminaLongClippingSeq extends IlluminaClippingSeq
	{
		IlluminaLongClippingSeq(String seq)
		{
			logger.infoln("Using Long Clipping Sequence: '" + seq + "'");

			this.seq = seq;

			long fullPack[] = packSeqInternal(seq, false);

			pack = new long[(fullPack.length + INTERLEAVE - 1) / INTERLEAVE];

			for (int i = 0; i < fullPack.length; i += INTERLEAVE)
				pack[i / INTERLEAVE] = fullPack[i];
		}

		public Integer readsSeqCompare(FastqRecord rec)
		{
			int seedMax = seedMaxMiss * 2;

			String recSequence = rec.getSequence();
			String clipSequence = seq;

			Set<Integer> offsetSet = new TreeSet<Integer>();

			long packRec[] = packSeqExternal(rec.getSequence());
			long packClip[] = getPack();

			int packRecMax = packRec.length - minSequenceOverlap;
			int packClipMax = packClip.length;

			for (int i = 0; i < packRecMax; i++)
				{
				long comboMask = calcSingleMask(packRec.length - i);

				for (int j = 0; j < packClipMax; j++)
					{
					int diff = Long.bitCount((packRec[i] ^ packClip[j]) & comboMask);

					if (diff <= seedMax)
						{
						int offset = i - j * INTERLEAVE;
						offsetSet.add(offset);
						}
					}
				}

			for (Integer offset : offsetSet)
				{
				int recCompLength = offset > 0 ? recSequence.length() - offset : recSequence.length();
				int clipCompLength = offset < 0 ? clipSequence.length() + offset : clipSequence.length();

				int compLength = recCompLength < clipCompLength ? recCompLength : clipCompLength;

				if (compLength > minSequenceOverlap)
					{
					float seqLikelihood = calculateDifferenceQuality(rec, clipSequence, compLength, offset);

					if (seqLikelihood >= minSequenceLikelihood)
						return offset;
					}
				}

			return null;
		}
	}

	private final static int BASE_A = 0x1;
	private final static int BASE_C = 0x4;
	private final static int BASE_G = 0x8;
	private final static int BASE_T = 0x2;

	public static long[] packSeqExternal(String seq)
	{
		long out[] = null;

		out = new long[seq.length()];

		long pack = 0;

		int offset = 0;

		for (int i = 0; i < 15; i++)
			{
			int tmp = 0;

			if (offset < seq.length())
				tmp = packCh(seq.charAt(offset), false);

			pack = (pack << 4) | tmp;
			offset++;
			}

		for (int i = 0; i < seq.length(); i++)
			{
			int tmp = 0;

			if (offset < seq.length())
				tmp = packCh(seq.charAt(offset), false);

			pack = (pack << 4) | tmp;
			out[i] = pack;

			offset++;
			}

		return out;
	}

	public static long calcSingleMask(int length)
	{
		/*
		 * if(length<16) return SHORT_MASKS[length];
		 * 
		 * return 0xFFFFFFFFFFFFFFFFL;
		 */
		long mask = 0xFFFFFFFFFFFFFFFFL;

		if (length < 16)
			mask <<= (16 - length) * 4;

		return mask;
	}

	public static long[] packSeqInternal(String seq, boolean reverse)
	{
		long out[] = null;

		if (!reverse)
			{
			out = new long[seq.length() - 15];

			long pack = 0;

			for (int i = 0; i < seq.length(); i++)
				{
				int tmp = packCh(seq.charAt(i), false);

				pack = (pack << 4) | tmp;

				if (i >= 15)
					out[i - 15] = pack;
				}
			}
		else
			{
			out = new long[seq.length() - 15];

			long pack = 0;

			for (int i = 0; i < seq.length(); i++)
				{
				long tmp = packCh(seq.charAt(i), true);

				pack = (pack >>> 4) | tmp << 60;

				if (i >= 15)
					out[i - 15] = pack;
				}

			}

		return out;
	}

	public static String unpack(long pack)
	{
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < 16; i++)
			{
			int tmp = (int) ((pack >>> 60) & 0xF);

			switch (tmp)
				{
				case BASE_A:
					sb.append("A");
					break;

				case BASE_C:
					sb.append("C");
					break;

				case BASE_G:
					sb.append("G");
					break;

				case BASE_T:
					sb.append("T");
					break;

				default:
					sb.append("[" + tmp + "]");
				}
			pack <<= 4;
			}

		return sb.toString();
	}

	private static int packCh(char ch, boolean rev)
	{
		if (!rev)
			{
			switch (ch)
				{
				case 'A':
					return BASE_A;
				case 'C':
					return BASE_C;
				case 'G':
					return BASE_G;
				case 'T':
					return BASE_T;
				}
			}
		else
			{
			switch (ch)
				{
				case 'A':
					return BASE_T;
				case 'C':
					return BASE_G;
				case 'G':
					return BASE_C;
				case 'T':
					return BASE_A;
				}
			}
		return 0;
	}
}
