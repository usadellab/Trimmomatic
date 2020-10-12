package org.usadellab.trimmomatic.fastq;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;
import org.itadaki.bzip2.BZip2InputStream;

public class FastqParser {

    private int phredOffset;
    private BufferedReader reader;
    private FastqRecord current;
    private AtomicInteger progress;
    private long fileLength;
    private long bytesRead;
    private EOL_TYPE eoltype;

    public FastqParser(int phredOffset) {
        this.phredOffset = phredOffset;
        this.progress = new AtomicInteger(0);
        this.bytesRead = 0;
    }

    public void parseOne() throws IOException {
        current = null;

        String name;
        String sequence;
        String comment;
        String quality;

        String line;

        line = reader.readLine();
        if (line == null) {
            progress.set(100);
            return;
        }
        bytesRead += line.length() + eoltype.getByteSize(); 

        if (line.startsWith("@")) {
            name = line.substring(1);
        } else {
            throw new RuntimeException("Invalid FASTQ name line: " + line);
        }

        sequence = reader.readLine();
        bytesRead += line.length() + eoltype.getByteSize(); 

        line = reader.readLine();
        bytesRead += line.length() + eoltype.getByteSize(); 
        if (line.startsWith("+")) {
            comment = line.substring(1);
        } else {
            throw new RuntimeException("Invalid FASTQ comment line: " + line);
        }

        quality = reader.readLine();
        bytesRead += line.length() + eoltype.getByteSize(); 

        current = new FastqRecord(name, sequence, comment, quality, phredOffset);
        progress.set((int) (((float) bytesRead / fileLength) * 100));
    }

    public int getProgress() {
        return this.progress.get();
    }

    public void parse(File file) throws IOException {
        String name = file.getName();
        fileLength = file.length();
        eoltype = guessEOLType(file);


        if (name.toLowerCase().endsWith(".gz")) {
            reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new BufferedInputStream(new FileInputStream(file), 1000000))));
        } else if (name.toLowerCase().endsWith(".bz2")) {
            reader = new BufferedReader(new InputStreamReader(new BZip2InputStream(new BufferedInputStream(new FileInputStream(file), 1000000), false)));
        } else if (name.toLowerCase().endsWith(".zip")) {
            reader = new BufferedReader(new InputStreamReader(new ZipInputStream(new BufferedInputStream(new FileInputStream(file), 1000000))));
        } else {
            reader = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(file), 1000000)));
        }

        parseOne();
    }

    public void close() throws IOException {
        reader.close();
    }

    public boolean hasNext() {
        return current != null;
    }

    public FastqRecord next() throws IOException {
        FastqRecord current = this.current;
        parseOne();

        return current;
    }

    public static enum EOL_TYPE {

        WINDOWS(2),
        UNIX(1),
        MACINTOSH(1),
        UNKNOWN(0);
        private int byteSize;

        private EOL_TYPE(int size) {
            this.byteSize = size;
        }

        public int getByteSize() {
            return this.byteSize;
        }
    }

    public static EOL_TYPE guessEOLType(File f) throws FileNotFoundException, IOException {
        InputStreamReader isr = new InputStreamReader(new FileInputStream(f));
        char symbol;
        while (isr.ready()) {
            symbol = (char) isr.read();
            if (symbol == '\n') {
                isr.close();
                return EOL_TYPE.UNIX;
            }
            if (symbol == '\r') {
                if (isr.read() == '\n') {
                    isr.close();
                    return EOL_TYPE.WINDOWS;
                } else {
                    isr.close();
                    return EOL_TYPE.MACINTOSH;
                }
            }
        }
        return EOL_TYPE.UNKNOWN;
    }
}
