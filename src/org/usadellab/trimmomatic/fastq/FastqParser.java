package org.usadellab.trimmomatic.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.ZipInputStream;

import org.itadaki.bzip2.BZip2InputStream;
import org.usadellab.trimmomatic.util.ConcatGZIPInputStream;
import org.usadellab.trimmomatic.util.PositionTrackingInputStream;

public class FastqParser {

    private int phredOffset;
    private PositionTrackingInputStream posTrackInputStream;
    private BufferedReader reader;
    private FastqRecord current;
    private long fileLength;

    private AtomicBoolean atEOF;
    
    public FastqParser(int phredOffset) {
        this.phredOffset = phredOffset;
        this.atEOF=new AtomicBoolean();
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
        	atEOF.set(true);
            return;
        }
        
        if (line.charAt(0)=='@') {
            name = line.substring(1);
        } else {
            throw new RuntimeException("Invalid FASTQ name line: " + line);
        }

        sequence = reader.readLine();

        line = reader.readLine();

        if (line.charAt(0)=='+') {
            comment = line.substring(1);
        } else {
            throw new RuntimeException("Invalid FASTQ comment line: " + line);
        }

        quality = reader.readLine();        
        current = new FastqRecord(name, sequence, comment, quality, phredOffset);
    }

    public int getProgress() {
    	if(atEOF.get())
    		return 100;
    	
    	long bytesRead=posTrackInputStream.getPosition();
    	
    	return (int)(((float) bytesRead / fileLength) * 100);    
    }

    public void parse(File file) throws IOException {
        String name = file.getName();
        fileLength = file.length();
        
        posTrackInputStream=new PositionTrackingInputStream(new FileInputStream(file));
        
        InputStream contentInputStream=posTrackInputStream;
        
        if (name.toLowerCase().endsWith(".gz")) {
            contentInputStream=new ConcatGZIPInputStream(posTrackInputStream);
        } else if (name.toLowerCase().endsWith(".bz2")) {
            contentInputStream=new BZip2InputStream(posTrackInputStream, false);
        } else if (name.toLowerCase().endsWith(".zip")) {
            contentInputStream=new ZipInputStream(posTrackInputStream);
        }
        
        reader=new BufferedReader(new InputStreamReader(contentInputStream), 32768);
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

}
