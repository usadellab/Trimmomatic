package org.usadellab.trimmomatic.trim;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * UMIINDEXMERGE:<indexFastqFile>[:<separator>]
 *
 * Merges a UMI carried in a separate index read (I2) FASTQ into the read
 * name, for Illumina UMI-adapter designs where the UMI physically replaces
 * the i7 index and is sequenced as its own index read rather than being
 * embedded in R1/R2 (e.g. IDT xGen UDI-UMI adapters, see
 * https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-library-preparation/ngs-adapters-indexing-primers/adapters-indexing-primers-for-illumina).
 * Every other UMI*-family trimmer assumes the UMI is embedded in the
 * sequence being trimmed; this one is for the case where it structurally
 * isn't, so a full 3rd FASTQ file has to be joined in by read ID.
 *
 * <indexFastqFile> is read completely into memory at construction (plain
 * text or .gz/.bz2/.zip, same as UMIDROPLETCORRECT/UMIDIMERCORRECT's
 * whitelist loading) into a read-ID -> UMI-sequence map. This is simpler and
 * safer than threading a synchronised second stream through Trimmomatic's
 * pipeline: Trimmomatic processes read blocks across multiple worker threads
 * concurrently, so "read the next line of the index file" per record would
 * race across threads and silently mismatch pairs. A read-ID keyed map has
 * no such ordering dependency, so correctness doesn't rely on block/thread
 * scheduling at all.
 *
 * The read ID used as the join key is the read name up to (but not
 * including) the first whitespace, matching how Illumina read names are
 * shared verbatim across R1/R2/I2 for the same template (only the space-
 * delimited description after it differs, e.g. "1:N:0:..." vs "2:N:0:...").
 *
 * Unlike UMIEXTRACT/UMISPLIT/UMIDROPLETCORRECT/UMIDIMERCORRECT/
 * UMIRHAPSODYCORRECT, this trimmer is SAFE in symmetric paired-end mode and
 * needs no processRecords() override: R1 and R2 of the same template share
 * the same read ID and therefore look up the exact same UMI from the index
 * map, so tagging each mate independently can never desynchronise their
 * names the way extracting a UMI from each mate's own leading bases would.
 *
 * A read whose ID has no entry in the index map is dropped (the index read
 * for that template is missing or the files are out of sync).
 *
 * The read is renamed as: @original_name<separator>UMI:<index read bases>,
 * matching UMIEXTRACT's own "UMI:"-labelled tagging convention, since this
 * is also a bare extract with no whitelist to correct against.
 *
 * Example:
 *   UMIINDEXMERGE:sample_I2.fastq.gz
 *
 * Path-vs-separator ambiguity: unlike UMIDROPLETCORRECT/UMIDIMERCORRECT,
 * there is no trailing numeric field here to tell a real trailing separator
 * apart from a Windows drive-letter colon ("C:\...") splitting the path
 * itself into an extra token. "Is the last token an integer" doesn't apply
 * when there's no integer field at all, so it's resolved instead by checking whether
 * the last colon-split token contains a path separator character ('/' or
 * '\\'): a real trailing <separator> argument realistically never contains
 * one, while any token that is actually part of a path very likely does
 * (either the path itself, or the token before an extensionless final
 * component). If the last token has no slash, it's treated as the
 * <separator>; otherwise every token is rejoined as the path and the
 * separator defaults to "_".
 */
public class UmiIndexMergeTrimmer extends AbstractSingleRecordTrimmer {
    private final String separator;
    private final Map<String, String> umiByReadId;

    public UmiIndexMergeTrimmer(String args) throws IOException {
        String[] arg = args.split(":");
        if (arg.length < 1 || arg[0].isEmpty())
            throw new IllegalArgumentException("UMIINDEXMERGE requires <indexFastqFile>, got: " + args);

        // See class javadoc "Path-vs-separator ambiguity" for why this check
        // (last token contains no path-separator character) is used instead of
        // the "is the last token an integer" heuristic the other UMI*Correct
        // trimmers use. There is no numeric trailing field here to key off.
        String lastToken = arg[arg.length - 1];
        boolean hasSeparator = arg.length > 1 && lastToken.indexOf('/') < 0 && lastToken.indexOf('\\') < 0;
        String path = hasSeparator
                ? String.join(":", java.util.Arrays.copyOfRange(arg, 0, arg.length - 1))
                : String.join(":", arg);
        separator = hasSeparator ? lastToken : "_";

        umiByReadId = loadIndexReads(path);
    }

    private static Map<String, String> loadIndexReads(String path) throws IOException {
        Map<String, String> map = new HashMap<>();
        FastqParser parser = new FastqParser(33);
        try {
            parser.open(new File(path));
            while (parser.hasNext()) {
                FastqRecord rec = parser.next();
                if (rec == null)
                    continue;
                map.put(readId(rec.getName()), rec.getSequence());
            }
        } finally {
            parser.close();
        }

        if (map.isEmpty())
            throw new IllegalArgumentException("UMIINDEXMERGE index FASTQ file is empty or unreadable: " + path);

        return map;
    }

    private static String readId(String name) {
        int space = name.indexOf(' ');
        return space < 0 ? name : name.substring(0, space);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        String umi = umiByReadId.get(readId(in.getName()));
        if (umi == null)
            return null;

        String newName = in.getName() + separator + "UMI:" + umi;
        return new FastqRecord(in, 0, in.getLength(), newName);
    }
}
