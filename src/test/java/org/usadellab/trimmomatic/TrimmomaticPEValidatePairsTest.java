package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;

import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;

/**
 * -validatePairs must stop the run when the mate names do not match, so a mispaired input
 * cannot finish with exit code 0.
 */
public class TrimmomaticPEValidatePairsTest {

	@TempDir
	Path tempDir;

	private File writeFastq(String filename, String... names) throws Exception {
		File f = tempDir.resolve(filename).toFile();
		try (FileWriter fw = new FileWriter(f)) {
			for (String name : names)
				fw.write("@" + name + "\nACGTACGTAC\n+\nIIIIIIIIII\n");
		}
		return f;
	}

	private void run(File r1, File r2, int threads) throws Exception {
		TrimmomaticPE tm = new TrimmomaticPE(new Logger(false, false, false));
		tm.process(r1, r2, false, tempDir.resolve("1P.fq").toFile(), tempDir.resolve("1U.fq").toFile(),
				tempDir.resolve("2P.fq").toFile(), tempDir.resolve("2U.fq").toFile(), new Trimmer[0], 33, null, null,
				true, null, null, threads, false, 0);
	}

	@ParameterizedTest
	@ValueSource(ints = { 1, 4 })
	@Timeout(30)
	public void testMismatchedPairsFail(int threads) throws Exception {
		File r1 = writeFastq("r1.fq", "read1/1", "read2/1", "read3/1");
		File r2 = writeFastq("r2.fq", "read1/2", "read3/2", "read2/2");

		assertThrows(RuntimeException.class, () -> run(r1, r2, threads));
	}

	@ParameterizedTest
	@ValueSource(ints = { 1, 4 })
	@Timeout(30)
	public void testMissingMateFails(int threads) throws Exception {
		File r1 = writeFastq("r1.fq", "read1/1", "read2/1", "read3/1");
		File r2 = writeFastq("r2.fq", "read1/2", "read2/2");

		assertThrows(RuntimeException.class, () -> run(r1, r2, threads));
	}

	@ParameterizedTest
	@ValueSource(ints = { 1, 4 })
	@Timeout(30)
	public void testMatchingPairsPass(int threads) throws Exception {
		File r1 = writeFastq("r1.fq", "read1/1", "read2/1", "read3/1");
		File r2 = writeFastq("r2.fq", "read1/2", "read2/2", "read3/2");

		assertDoesNotThrow(() -> run(r1, r2, threads));
	}
}
