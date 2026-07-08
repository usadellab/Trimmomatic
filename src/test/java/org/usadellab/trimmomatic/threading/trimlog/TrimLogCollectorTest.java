package org.usadellab.trimmomatic.threading.trimlog;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class TrimLogCollectorTest {

	@Test
	public void testParasiteCollectorWritesExpectedContent() throws Exception {
		File trimLog = File.createTempFile("trimlog", ".log");
		trimLog.deleteOnExit();

		TrimLogCollector collector = TrimLogCollector.makeTrimLogCollector(false, 1, trimLog, new ExceptionHolder());

		BlockOfRecords block = new BlockOfRecords(null, null);
		block.setTrimLogRecs(List.of(new TrimLogRecord("seq1", 150, 0, 150, 0),
				new TrimLogRecord("seq2", 100, 5, 95, 10)));

		collector.put(block);
		collector.close();

		String content = Files.readString(trimLog.toPath());
		String expected = "seq1 150 0 150 0\n" + "seq2 100 5 95 10\n";

		assertEquals(expected, content);
	}

	@Test
	public void testCollectorFlushesManyRecordsOnClose() throws Exception {
		File trimLog = File.createTempFile("trimlog", ".log");
		trimLog.deleteOnExit();

		TrimLogCollector collector = TrimLogCollector.makeTrimLogCollector(false, 1, trimLog, new ExceptionHolder());

		int recordCount = 10_000;
		TrimLogRecord[] recs = new TrimLogRecord[recordCount];
		for (int i = 0; i < recordCount; i++)
			recs[i] = new TrimLogRecord("seq" + i, 150, 0, 150, 0);

		BlockOfRecords block = new BlockOfRecords(null, null);
		block.setTrimLogRecs(List.of(recs));

		collector.put(block);
		collector.close();

		long lineCount;
		try (var lines = Files.lines(trimLog.toPath())) {
			lineCount = lines.count();
		}

		assertEquals(recordCount, lineCount);
	}
}
