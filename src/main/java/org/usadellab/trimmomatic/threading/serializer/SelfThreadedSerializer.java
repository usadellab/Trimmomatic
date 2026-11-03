package org.usadellab.trimmomatic.threading.serializer;

import java.io.IOException;
import java.io.OutputStream;

import org.usadellab.trimmomatic.threading.ExceptionHolder;
import org.usadellab.trimmomatic.util.compression.BlockData;

public class SelfThreadedSerializer extends Serializer implements Runnable {
	protected ExceptionHolder exceptionHolder;

	private Thread thread;

	public SelfThreadedSerializer(OutputStream stream, SerializedBlockQueue sbq, ExceptionHolder exceptionHolder)
			throws IOException {
		super(stream, sbq);

		this.exceptionHolder = exceptionHolder;

		thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}

	public void pollWritable() throws Exception {
		exceptionHolder.rethrow();
	}

	public void close() throws Exception {
		try {
			while (thread.isAlive()) {
				exceptionHolder.rethrow();
				thread.join(100);
			}
		} catch (Exception e) {
			// The run failed, the output is incomplete anyway. Release the file handle
			// without writing a trailer.
			try {
				stream.close();
			} catch (Exception closeFailure) {
				e.addSuppressed(closeFailure);
			}
			throw e;
		}

		super.close();
	}

	@Override
	public void run() {
		try {
			BlockData blockData = sbq.takeOutput().get();

			while (blockData != null) {
				stream.writeBlock(blockData);
				blockData = sbq.takeOutput().get();
			}
		} catch (Exception e) {
			Exception pe = new Exception("Serializer Exception", e);
			exceptionHolder.setException(pe);
		} finally {
			setCompleted();
		}
	}

}
