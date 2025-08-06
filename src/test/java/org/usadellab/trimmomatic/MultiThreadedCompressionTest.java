package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.doAnswer;
import static org.mockito.Mockito.mock;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ExceptionHolder;
import org.usadellab.trimmomatic.threading.pipeline.Pipeline;
import org.usadellab.trimmomatic.threading.serializer.Serializer;
import org.usadellab.trimmomatic.trim.Trimmer;

public class MultiThreadedCompressionTest {

    @Test
    void testThreadedPipelineExecution() throws Exception {
        // 1. ARRANGE: Set up the multi-threaded pipeline and a block of work.

        // A latch to allow the main thread to wait for the background thread to finish.
        final CountDownLatch latch = new CountDownLatch(1);

        ExceptionHolder exceptionHolder = new ExceptionHolder();
        
        // Create a pipeline with 2 threads.
        Pipeline pipeline = Pipeline.makePipeline(2, exceptionHolder);

        // Create mock objects for dependencies that aren't the focus of this test.
        Trimmer[] trimmers = new Trimmer[0]; // No actual trimming needed
        
        // Mock the Serializer to signal when it's done.
        Serializer serializer = mock(Serializer.class);
        doAnswer(invocation -> {
            latch.countDown(); // Signal that the work is complete.
            return null;
        }).when(serializer).pollCompressible();
        
        List<Serializer> serializers = new ArrayList<>();
        serializers.add(serializer);

        // Create some initial data to be processed.
        List<FastqRecord> records = new ArrayList<>();
        records.add(new FastqRecord("read1", "GATTACA", "comment1", "!!!!!!!", 33));
        BlockOfRecords bor = new BlockOfRecords(records, null);

        // Create the 'BlockOfWork' which is the task to be run on a separate thread.
        BlockOfWork work = new BlockOfWork(null, trimmers, bor, true, false, false, serializers, exceptionHolder);
        
        
        // 2. ACT: Submit the work to the pipeline.
        
        Future<BlockOfRecords> future = pipeline.submit(work);

        // 3. ASSERT: Verify that the task completed successfully.

        // Wait for the future to complete. This ensures the 'call()' method finished.
        BlockOfRecords result = future.get(5, TimeUnit.SECONDS);
        assertNotNull(result, "The future should return a valid BlockOfRecords object.");

        // Wait for the serializer to be called. This ensures the full 'process()' method ran.
        boolean completedInTime = latch.await(5, TimeUnit.SECONDS);
        assertTrue(completedInTime, "The background task should complete and call the serializer.");
        
        // Check if any exceptions were caught by the pipeline.
        exceptionHolder.rethrow();
        
        // Clean up the pipeline.
        pipeline.close();
    }
}