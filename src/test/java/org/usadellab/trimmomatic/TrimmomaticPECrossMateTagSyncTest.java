package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.trim.HeadCropTrimmer;
import org.usadellab.trimmomatic.trim.MinLenTrimmer;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.trim.UmiExtractTrimmer;
import org.usadellab.trimmomatic.trim.UmiSplitTrimmer;
import org.usadellab.trimmomatic.util.Logger;

/**
 * Real-world single-cell tooling (umi_tools extract, and anything that reads
 * a UMI/cell-barcode back out of an aligned BAM's QNAME) writes the extracted
 * tag onto BOTH mates' names, not just the mate it was extracted from --
 * because the technical mate (e.g. 10x R1) is normally discarded before
 * alignment, only the biological mate (R2) survives into a BAM, and a tag
 * that lives only on the discarded mate's name is invisible to anything
 * downstream that reads it back out post-alignment.
 *
 * BlockOfWork's per-mate branch now mirrors any appended name tag onto the
 * OTHER mate automatically -- detected generically (any name-tagging step
 * only ever appends, never replaces), so this works for UMIEXTRACT, UMISPLIT
 * and BARCODECORRECT without any of those classes needing to know about the
 * other mate at all.
 */
public class TrimmomaticPECrossMateTagSyncTest {

    @TempDir
    Path tempDir;

    private static String fqRecord(String name, String seq) {
        String qual = "I".repeat(seq.length());
        return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
    }

    private File writeFastq(String filename, String... records) throws Exception {
        File f = tempDir.resolve(filename).toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String r : records) fw.write(r);
        }
        return f;
    }

    private List<String> readNames(File f) throws Exception {
        List<String> names = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(f))) {
            String line;
            int lineNo = 0;
            while ((line = br.readLine()) != null) {
                if (lineNo % 4 == 0) names.add(line.substring(1)); // strip leading '@'
                lineNo++;
            }
        }
        return names;
    }

    private void run(File r1in, File r2in, File r1p, File r1u, File r2p, File r2u,
                     Trimmer[] trimmers1, Trimmer[] trimmers2) throws Exception {
        Logger logger = new Logger(false, false, false);
        TrimmomaticPE tm = new TrimmomaticPE(logger);
        tm.process(r1in, r2in, false, r1p, r1u, r2p, r2u, new Trimmer[0], trimmers1, trimmers2,
                33, null, null, false, null, null, 1, false, 0);
    }

    @Test
    public void testUmiSplitOnMate1_tagMirroredOntoMate2() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "AAAACCCC" + "TTTT")); // 4bp CB + 4bp UMI + payload
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "GGGGGGGGGGGGGGGG"));
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u,
                new Trimmer[] { new UmiSplitTrimmer("4:4") },
                new Trimmer[0]);

        List<String> r1Names = readNames(r1p);
        List<String> r2Names = readNames(r2p);
        assertEquals(1, r1Names.size());
        assertEquals(1, r2Names.size());

        assertEquals("p1/1_CB:AAAA_UMI:CCCC", r1Names.get(0), "mate 1 keeps its own tagged name");
        assertEquals("p1/2_CB:AAAA_UMI:CCCC", r2Names.get(0),
                "mate 2 must have the SAME tag mirrored onto its own (otherwise-untouched) name");
    }

    @Test
    public void testUmiExtractOnMate1_tagMirroredOntoMate2() throws Exception {
        // Generalisation check: this mechanism isn't UMISPLIT/BARCODECORRECT-specific --
        // any name-appending step gets mirrored, including the pre-existing UMIEXTRACT.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGT" + "TTTT")); // 8bp UMI + payload
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "GGGGGGGGGGGGGGGG"));
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u,
                new Trimmer[] { new UmiExtractTrimmer("8") },
                new Trimmer[0]);

        List<String> r1Names = readNames(r1p);
        List<String> r2Names = readNames(r2p);

        assertEquals("p1/1_UMI:ACGTACGT", r1Names.get(0));
        assertEquals("p1/2_UMI:ACGTACGT", r2Names.get(0),
                "UMIEXTRACT's tag must also mirror onto mate 2, same mechanism as UMISPLIT");
    }

    @Test
    public void testNeitherMateRenames_namesUntouched() throws Exception {
        // Regression: the common case (quality trim only, no renaming step anywhere)
        // must be completely unaffected -- this is the exact scenario Phase 0's own
        // tests already covered, confirming this change doesn't disturb it.
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "ACGTACGTACGT"));
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "TTTTGGGGCCCC"));
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u,
                new Trimmer[] { new HeadCropTrimmer("4") },
                new Trimmer[] { new MinLenTrimmer(1) });

        assertEquals("p1/1", readNames(r1p).get(0), "no tagging step used -- name must be completely untouched");
        assertEquals("p1/2", readNames(r2p).get(0), "no tagging step used -- name must be completely untouched");
    }

    @Test
    public void testBothMatesTagIndependently_unionAppliedToBoth() throws Exception {
        File r1in = writeFastq("r1.fastq", fqRecord("p1/1", "AAAA" + "TTTT")); // 4bp UMI + payload
        File r2in = writeFastq("r2.fastq", fqRecord("p1/2", "CCCC" + "GGGG")); // 4bp UMI + payload
        File r1p = tempDir.resolve("r1p.fastq").toFile();
        File r1u = tempDir.resolve("r1u.fastq").toFile();
        File r2p = tempDir.resolve("r2p.fastq").toFile();
        File r2u = tempDir.resolve("r2u.fastq").toFile();

        run(r1in, r2in, r1p, r1u, r2p, r2u,
                new Trimmer[] { new UmiExtractTrimmer("4") },
                new Trimmer[] { new UmiExtractTrimmer("4") });

        assertEquals("p1/1_UMI:AAAA_UMI:CCCC", readNames(r1p).get(0),
                "mate 1 keeps its own name, gains BOTH tags (its own plus mate 2's)");
        assertEquals("p1/2_UMI:AAAA_UMI:CCCC", readNames(r2p).get(0),
                "mate 2 keeps its own name, gains BOTH tags (its own plus mate 1's) -- identical suffix to mate 1's");
    }
}
