package org.usadellab.trimmomatic;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fasta.FastaRecord;

public class TrimmomaticTest {

	@Test
	void testReverseComplement() {
		// 1. Arrange: Define the input and expected output
		String originalSequence = "AGTCG";
		String expectedReverseComplement = "CGACT";

		// 2. Act: Call the method you want to test
		String actualReverseComplement = FastaRecord.getComplementSequence(originalSequence);

		// 3. Assert: Check if the actual result matches the expected result
		assertEquals(expectedReverseComplement, actualReverseComplement,
				"The reverse complement was not calculated correctly.");
	}
}