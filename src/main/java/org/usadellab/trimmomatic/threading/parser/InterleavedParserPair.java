package org.usadellab.trimmomatic.threading.parser;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

/**
 * Reads a single interleaved FASTQ source (R1 and R2 alternating) and
 * exposes two Parser views — one for the odd-indexed records (R1) and
 * one for the even-indexed records (R2) — that can be plugged directly
 * into the standard TrimmomaticPE pipeline.
 *
 * A single background thread drives all I/O; both half-parsers share its
 * output via independent blocking queues, so the main processing loop
 * never blocks on I/O.
 */
public class InterleavedParserPair implements AutoCloseable {

    private final ArrayBlockingQueue<List<FastqRecord>> r1Queue;
    private final ArrayBlockingQueue<List<FastqRecord>> r2Queue;
    private final ExceptionHolder exceptionHolder;
    private final Thread readerThread;
    private final FastqParser source;

    private final HalfParser r1Parser;
    private final HalfParser r2Parser;

    public InterleavedParserPair(FastqParser source, int buffers,
            ExceptionHolder exceptionHolder) {
        this.source = source;
        this.exceptionHolder = exceptionHolder;
        this.r1Queue = new ArrayBlockingQueue<>(buffers);
        this.r2Queue = new ArrayBlockingQueue<>(buffers);
        this.r1Parser = new HalfParser(r1Queue, exceptionHolder);
        this.r2Parser = new HalfParser(r2Queue, exceptionHolder);

        readerThread = new Thread(this::readAll, "InterleavedParser");
        readerThread.setDaemon(true);
        readerThread.start();
    }

    public Parser getR1Parser() { return r1Parser; }
    public Parser getR2Parser() { return r2Parser; }

    public void close() throws Exception {
        while (readerThread.isAlive()) {
            exceptionHolder.rethrow();
            readerThread.join(100);
        }
        // Check once more after the thread has died — the exception might have
        // been stored after the last isAlive() check in the loop.
        try {
            exceptionHolder.rethrow();
        } finally {
            // Always release the file handle, even when rethrowing.
            source.close();
        }
    }

    // -----------------------------------------------------------------------

    private void readAll() {
        try {
            while (source.hasNext()) {
                // Cap initial list capacity for long-read workloads.
                List<FastqRecord> r1 = new ArrayList<>(Math.min(Parser.BLOCK_MAX_RECORDS, 1024));
                List<FastqRecord> r2 = new ArrayList<>(Math.min(Parser.BLOCK_MAX_RECORDS, 1024));

                // Read up to BLOCK_MAX_RECORDS pairs OR BLOCK_MAX_BYTES total bytes.
                int  pairs      = 0;
                long blockBytes = 0;
                while (source.hasNext()
                        && pairs      <  Parser.BLOCK_MAX_RECORDS
                        && blockBytes <  Parser.BLOCK_MAX_BYTES) {
                    FastqRecord rec1 = source.next();
                    if (!source.hasNext())
                        throw new RuntimeException(
                                "Interleaved FASTQ has an odd number of records — " +
                                "R1 record '" + rec1.getName() + "' has no R2 partner.");
                    FastqRecord rec2 = source.next();
                    r1.add(rec1);
                    r2.add(rec2);
                    blockBytes += rec1.getRecordLength() + rec2.getRecordLength();
                    pairs++;
                }

                r1Queue.put(r1);
                r2Queue.put(r2);
            }
        } catch (Exception e) {
            exceptionHolder.setException(new Exception("InterleavedParser Exception", e));
        } finally {
            // EOF sentinels
            try {
                r1Queue.put(new ArrayList<>());
                r2Queue.put(new ArrayList<>());
            } catch (InterruptedException e) {
                exceptionHolder.setException(e);
            }
        }
    }

    // -----------------------------------------------------------------------

    /**
     * A lightweight Parser that simply drains a pre-filled blocking queue.
     * All I/O is handled by the InterleavedParserPair reader thread; this
     * class just exposes the standard Parser.poll() interface.
     */
    private static final class HalfParser extends Parser {
        private final ArrayBlockingQueue<List<FastqRecord>> queue;
        private final ExceptionHolder exceptionHolder;
        private final AtomicBoolean closed = new AtomicBoolean(false);

        HalfParser(ArrayBlockingQueue<List<FastqRecord>> queue,
                ExceptionHolder exceptionHolder) {
            // Pass null — this parser never calls parseBlock() or super.close().
            super(null);
            this.queue = queue;
            this.exceptionHolder = exceptionHolder;
        }

        @Override
        public List<FastqRecord> poll() throws Exception {
            exceptionHolder.rethrow();
            return queue.poll(100, TimeUnit.MILLISECONDS);
        }

        @Override
        public void close() {
            closed.set(true);
            // Actual source close is handled by InterleavedParserPair.close().
        }
    }
}
