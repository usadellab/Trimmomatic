package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;
import org.usadellab.trimmomatic.util.compression.BlockData;
import org.usadellab.trimmomatic.util.compression.BlockOutputStream;
import org.usadellab.trimmomatic.util.compression.CompressionFormat;
import org.usadellab.trimmomatic.util.compression.ParallelCompressor;
import org.usadellab.trimmomatic.util.compression.UncompressedBlockData;

public class CompressionTest {

    // Updated to also provide 'null' for the default compression level
    static Stream<Arguments> compressorProvider() {
        List<String> formats = Arrays.asList("test.gz", "test.bz2");
        List<Integer> levels = Arrays.asList(1, 5, 9); // Specific levels to test
        
        List<Arguments> testCases = new ArrayList<>();
        for (String format : formats) {
            // Add test cases for specific levels
            for (Integer level : levels) {
                testCases.add(Arguments.of(format, level));
            }
            // Add test case for the default (null) level
            testCases.add(Arguments.of(format, null));
        }
        
        return testCases.stream();
    }

    // Updated test signature to use Integer to allow for null values
    @ParameterizedTest
    @MethodSource("compressorProvider")
    void testParallelCompression(String outputFileName, Integer compressLevel) throws Exception {
        // 1. ARRANGE: Set up test data and compressor
        ParallelCompressor compressor = CompressionFormat.parallelCompressorForSerializing(null, outputFileName, compressLevel);

        List<FastqRecord> records = new ArrayList<>();
        int phredOffset = 33;
        records.add(new FastqRecord("read1", "GATTACA", "comment1", "!!!!!!!", phredOffset));
        records.add(new FastqRecord("read2", "TACCAGA", "comment2", "FFFFFFF", phredOffset));

        UncompressedBlockData uncompressedBlock = new UncompressedBlockData(records);
        
        ByteArrayOutputStream compressedByteStream = new ByteArrayOutputStream();

        // 2. ACT: Run the full compression lifecycle
        BlockOutputStream blockOutStream = compressor.wrapAndWriteHeader(compressedByteStream);
        compressor.updateChecksumPreCompression(uncompressedBlock);
        BlockData compressedBlock = compressor.compress(null, uncompressedBlock);
        compressor.updateChecksumPostCompression(compressedBlock);
        blockOutStream.writeBlock(compressedBlock);
        compressor.writeTrailer(blockOutStream);
        blockOutStream.close();

        // 3. ASSERT: Decompress and verify the output
        byte[] compressedBytes = compressedByteStream.toByteArray();
        InputStream decompressorStream = CompressionFormat.wrapStreamForParsing(new ByteArrayInputStream(compressedBytes), outputFileName);
        
        ByteArrayOutputStream result = new ByteArrayOutputStream();
        byte[] buffer = new byte[1024];
        int length;
        while ((length = decompressorStream.read(buffer)) != -1) {
            result.write(buffer, 0, length);
        }
        String decompressedData = result.toString(StandardCharsets.UTF_8.name());

        String expectedData = "@read1\nGATTACA\n+comment1\n!!!!!!!\n@read2\nTACCAGA\n+comment2\nFFFFFFF\n";

        // Updated assertion message to handle the null case
        String levelDescription = (compressLevel == null) ? "default" : "level " + compressLevel;
        assertEquals(expectedData, decompressedData, 
            "Decompressed data should match original for " + outputFileName + " at " + levelDescription);
    }
}