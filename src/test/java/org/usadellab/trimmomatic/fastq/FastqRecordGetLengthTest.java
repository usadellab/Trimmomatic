package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.lang.reflect.Field;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.trim.CropTrimmer;
import org.usadellab.trimmomatic.trim.HeadCropTrimmer;
import org.usadellab.trimmomatic.trim.MaxLenTrimmer;
import org.usadellab.trimmomatic.trim.MinLenTrimmer;
import org.usadellab.trimmomatic.trim.TailCropTrimmer;

/**
 * Verifies that FastqRecord.getLength() returns the correct value for plain
 * records, view records, and nested views, and that calling it on a view does
 * NOT trigger lazy materialisation of the sequence string.
 */
public class FastqRecordGetLengthTest {

    // --- helpers ---

    private FastqRecord makeRecord(int len) {
        StringBuilder seq = new StringBuilder(len);
        StringBuilder qual = new StringBuilder(len);
        for (int i = 0; i < len; i++) {
            seq.append("ACGT".charAt(i % 4));
            qual.append('I'); // phred33 quality 40
        }
        return new FastqRecord("read", seq.toString(), "", qual.toString(), 33);
    }

    /** Returns the value of the private {@code sequence} field via reflection. */
    private String rawSequenceField(FastqRecord rec) throws Exception {
        Field f = FastqRecord.class.getDeclaredField("sequence");
        f.setAccessible(true);
        return (String) f.get(rec);
    }

    // --- getLength() correctness ---

    @Test
    public void testGetLengthPlainRecord() {
        FastqRecord rec = makeRecord(150);
        assertEquals(150, rec.getLength());
    }

    @Test
    public void testGetLengthViewRecord() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 50, 60); // bases 50..109
        assertEquals(60, view.getLength());
    }

    @Test
    public void testGetLengthNestedView() {
        FastqRecord base = makeRecord(150);
        FastqRecord view1 = new FastqRecord(base, 10, 100); // bases 10..109
        FastqRecord view2 = new FastqRecord(view1, 5, 40);  // bases 15..54
        assertEquals(40, view2.getLength());
    }

    @Test
    public void testGetLengthClampedView() {
        // Requesting more bases than available should clamp, not throw.
        FastqRecord base = makeRecord(100);
        FastqRecord view = new FastqRecord(base, 80, 50); // only 20 bases available
        assertEquals(20, view.getLength());
    }

    @Test
    public void testGetLengthConsistentWithGetSequence() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 30, 80);
        // Both before and after materialisation getLength() must equal getSequence().length()
        assertEquals(view.getSequence().length(), view.getLength());
    }

    @Test
    public void testGetLengthEmptyView() {
        FastqRecord base = makeRecord(50);
        FastqRecord view = new FastqRecord(base, 50, 0); // zero-length trim
        assertEquals(0, view.getLength());
    }

    // --- no premature materialisation ---

    @Test
    public void testGetLengthDoesNotMaterialiseView() throws Exception {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 20, 80);

        // Before any access the sequence field must be null (it's a view).
        assertNull(rawSequenceField(view),
                "sequence field should be null on a fresh view record");

        // getLength() must not trigger materialisation.
        int len = view.getLength();
        assertEquals(80, len);
        assertNull(rawSequenceField(view),
                "getLength() must not materialise the sequence substring");

        // Only after getSequence() should the field become non-null.
        view.getSequence();
        assertNotNull(rawSequenceField(view),
                "sequence field should be non-null after getSequence()");
    }

    @Test
    public void testGetLengthDoesNotMaterialiseNestedView() throws Exception {
        FastqRecord base = makeRecord(150);
        FastqRecord view1 = new FastqRecord(base, 10, 130);
        FastqRecord view2 = new FastqRecord(view1, 10, 50);

        assertNull(rawSequenceField(view2));
        assertEquals(50, view2.getLength());
        assertNull(rawSequenceField(view2),
                "getLength() on nested view must not materialise the sequence");
    }

    // --- trimmer integration: views fed to length-only trimmers ---

    @Test
    public void testMinLenTrimmerOnViewKeeps() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 0, 80); // 80 bp view

        MinLenTrimmer trimmer = new MinLenTrimmer(80);
        assertNotNull(trimmer.processRecord(view)); // exactly at threshold → keep
    }

    @Test
    public void testMinLenTrimmerOnViewDrops() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 0, 79); // 79 bp view

        MinLenTrimmer trimmer = new MinLenTrimmer(80);
        assertNull(trimmer.processRecord(view)); // one short → drop
    }

    @Test
    public void testMaxLenTrimmerOnViewKeeps() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 0, 50); // 50 bp view

        MaxLenTrimmer trimmer = new MaxLenTrimmer(50);
        assertNotNull(trimmer.processRecord(view));
    }

    @Test
    public void testMaxLenTrimmerOnViewDrops() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 0, 51); // 51 bp view

        MaxLenTrimmer trimmer = new MaxLenTrimmer(50);
        assertNull(trimmer.processRecord(view));
    }

    @Test
    public void testCropTrimmerOnView() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 20, 100); // 100 bp view

        CropTrimmer trimmer = new CropTrimmer(60);
        FastqRecord result = trimmer.processRecord(view);

        assertNotNull(result);
        assertEquals(60, result.getLength());
        assertEquals(60, result.getSequence().length());
        // Content should be the first 60 bases of the view (base positions 20..79)
        assertEquals(base.getSequence().substring(20, 80), result.getSequence());
    }

    @Test
    public void testCropTrimmerOnViewShorterThanCrop() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 20, 40); // 40 bp view

        CropTrimmer trimmer = new CropTrimmer(100); // crop larger than view → pass through
        FastqRecord result = trimmer.processRecord(view);

        assertNotNull(result);
        assertEquals(40, result.getLength());
    }

    @Test
    public void testHeadCropTrimmerOnView() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 10, 100); // 100 bp view

        HeadCropTrimmer trimmer = new HeadCropTrimmer("15");
        FastqRecord[] results = trimmer.processRecords(new FastqRecord[]{ view });

        assertNotNull(results[0]);
        assertEquals(85, results[0].getLength());
        // Content: base positions 25..109
        assertEquals(base.getSequence().substring(25, 110), results[0].getSequence());
    }

    @Test
    public void testTailCropTrimmerOnView() {
        FastqRecord base = makeRecord(150);
        FastqRecord view = new FastqRecord(base, 10, 100); // 100 bp view

        TailCropTrimmer trimmer = new TailCropTrimmer("15");
        FastqRecord[] results = trimmer.processRecords(new FastqRecord[]{ view });

        assertNotNull(results[0]);
        assertEquals(85, results[0].getLength());
        // Content: base positions 10..94
        assertEquals(base.getSequence().substring(10, 95), results[0].getSequence());
    }
}
