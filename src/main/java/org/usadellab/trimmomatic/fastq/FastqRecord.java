package org.usadellab.trimmomatic.fastq;

public class FastqRecord {
	private String name;
	private String sequence;
	private String comment;
	private String quality;
	private String barcodeLabel = null;

	private int phredOffset;
	private int headPos;

	// View fields — non-null only when this record is a trimmed view of another.
	// Avoids substring() allocations for intermediate trimming steps; sequence and
	// quality are materialised lazily on first access via getSequence()/getQuality().
	private String rawSequence;
	private String rawQuality;
	private int viewOffset;
	private int viewLength;

	// Quality integer arrays — computed lazily and invalidated when phredOffset changes.
	private int[] qualityCacheRaw;
	private int[] qualityCacheZeroNs;

	public FastqRecord(String name, String sequence, String comment, String quality, int phredOffset) {
		this.name = name;
		this.sequence = sequence;
		this.comment = comment;
		this.quality = quality;

		this.phredOffset = phredOffset;
		headPos = 0;

		if (sequence.length() != quality.length())
			throw new RuntimeException(
					"Sequence and quality length don't match: '" + sequence + "' vs '" + quality + "'");
	}

	public FastqRecord(FastqRecord base, int headPos, int length) {
		if (headPos < 0)
			throw new RuntimeException("Attempting invalid trim on " + base.name + " with length "
					+ base.getLength() + ": Wanted " + headPos + " to " + (headPos + length));

		int availableLength = base.getLength();
		if (headPos + length > availableLength)
			length = availableLength - headPos;

		if (length < 0)
			throw new StringIndexOutOfBoundsException(
					"begin " + headPos + ", end " + (headPos + length) + ", length " + availableLength);

		// Chain directly to the root raw strings so layered views never form.
		this.rawSequence = (base.rawSequence != null) ? base.rawSequence : base.sequence;
		this.rawQuality  = (base.rawQuality  != null) ? base.rawQuality  : base.quality;
		this.viewOffset  = (base.rawSequence != null) ? (base.viewOffset + headPos) : headPos;
		this.viewLength  = length;

		this.name        = base.name;
		this.comment     = base.comment;
		this.phredOffset = base.phredOffset;
		this.headPos     = base.headPos + headPos;
		this.barcodeLabel = base.barcodeLabel;
	}

	/**
	 * View constructor with an explicit name override. Used by UmiExtractTrimmer
	 * to append the extracted UMI to the read name without materialising the
	 * sequence or quality substrings.
	 */
	public FastqRecord(FastqRecord base, int headPos, int length, String nameOverride) {
		this(base, headPos, length);
		this.name = nameOverride;
	}

	/**
	 * View constructor with explicit name and comment overrides. Used by
	 * LongReadTrimmer so that the '+' quality header matches the renamed '@'
	 * sequence header, keeping output valid for strict FASTQ parsers.
	 */
	public FastqRecord(FastqRecord base, int headPos, int length, String nameOverride, String commentOverride) {
		this(base, headPos, length);
		this.name = nameOverride;
		this.comment = commentOverride;
	}

	public FastqRecord(FastqRecord base, String sequence, String quality, int phredOffset) {
		this.sequence = sequence;
		this.quality = quality;
		this.name = base.name;
		this.comment = base.comment;
		this.headPos = base.headPos;
		this.phredOffset = phredOffset;

		this.barcodeLabel = base.barcodeLabel;
	}

	public static FastqRecord make(String name, String sequence, int quality) {
		StringBuilder qualityBuilder = new StringBuilder();
		for (int i = 0; i < sequence.length(); i++)
			qualityBuilder.append((char) (33 + quality));

		return new FastqRecord(name, sequence, name, qualityBuilder.toString(), 33);
	}

	public static FastqRecord make(String name, String sequence) {
		return make(name, sequence, 40);
	}

	public int getLength() {
		return (rawSequence != null) ? viewLength : sequence.length();
	}

	public String getName() {
		return name;
	}

	public String getSequence() {
		if (sequence == null)
			sequence = rawSequence.substring(viewOffset, viewOffset + viewLength);
		return sequence;
	}

	public String getBarcodeLabel() {
		return barcodeLabel;
	}

	public void setBarcodeLabel(String barcodeLabel) {
		this.barcodeLabel = barcodeLabel;
	}

	public String getComment() {
		return comment;
	}

	public String getQuality() {
		if (quality == null)
			quality = rawQuality.substring(viewOffset, viewOffset + viewLength);
		return quality;
	}

	public int getPhredOffset() {
		return phredOffset;
	}

	void setPhredOffset(int phredOffset) {
		this.phredOffset = phredOffset;
		this.qualityCacheRaw = null;
		this.qualityCacheZeroNs = null;
	}

	public int getHeadPos() {
		return headPos;
	}

	public int[] getQualityAsInteger(boolean zeroNs) {
		if (zeroNs) {
			if (qualityCacheZeroNs != null)
				return qualityCacheZeroNs;
		} else {
			if (qualityCacheRaw != null)
				return qualityCacheRaw;
		}

		int[] arr;

		if (rawQuality != null) {
			arr = new int[viewLength];
			for (int i = 0; i < viewLength; i++) {
				if (zeroNs && rawSequence.charAt(viewOffset + i) == 'N')
					arr[i] = 0;
				else
					arr[i] = rawQuality.charAt(viewOffset + i) - phredOffset;
			}
		} else {
			arr = new int[quality.length()];
			for (int i = 0; i < quality.length(); i++) {
				if (zeroNs && sequence.charAt(i) == 'N')
					arr[i] = 0;
				else
					arr[i] = quality.charAt(i) - phredOffset;
			}
		}

		if (zeroNs)
			qualityCacheZeroNs = arr;
		else
			qualityCacheRaw = arr;

		return arr;
	}

	private static int RECORD_ADDED_LENGTH = 10;

	public int getRecordLength() {
		int seqLen = getLength();
		return this.name.length() + seqLen + this.comment.length() + seqLen + RECORD_ADDED_LENGTH;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((comment == null) ? 0 : comment.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + phredOffset;
		String q = getQuality();
		result = prime * result + ((q == null) ? 0 : q.hashCode());
		String s = getSequence();
		result = prime * result + ((s == null) ? 0 : s.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		FastqRecord other = (FastqRecord) obj;
		if (comment == null) {
			if (other.comment != null)
				return false;
		} else if (!comment.equals(other.comment))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (phredOffset != other.phredOffset)
			return false;
		String q = getQuality();
		String oq = other.getQuality();
		if (q == null) {
			if (oq != null)
				return false;
		} else if (!q.equals(oq))
			return false;
		String s = getSequence();
		String os = other.getSequence();
		if (s == null) {
			if (os != null)
				return false;
		} else if (!s.equals(os))
			return false;
		return true;
	}

}
