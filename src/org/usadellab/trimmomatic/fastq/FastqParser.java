package org.usadellab.trimmomatic.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayDeque;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.ZipInputStream;

import org.itadaki.bzip2.BZip2InputStream;
import org.usadellab.trimmomatic.util.ConcatGZIPInputStream;
import org.usadellab.trimmomatic.util.PositionTrackingInputStream;

public class FastqParser {

	private static final int PREREAD_COUNT=10000;

    private int phredOffset;
    private ArrayDeque<FastqRecord> deque;
    int qualHistogram[];
    int patternHistogram[];
    
    private PositionTrackingInputStream posTrackInputStream;
    private BufferedReader reader;
    private FastqRecord current;
    private long fileLength;

    private AtomicBoolean atEOF;
    
    public FastqParser(int phredOffset) {
        this.phredOffset = phredOffset;
        deque=new ArrayDeque<FastqRecord>(PREREAD_COUNT);
        
        this.atEOF=new AtomicBoolean();
    }

    public void setPhredOffset(int phredOffset)
    {
    	this.phredOffset=phredOffset;
    	
    	if(current!=null)
    		current.setPhredOffset(phredOffset);
    }
    
    public void parseOne() throws IOException 
    {
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
        if(sequence==null)
        	throw new RuntimeException("Missing sequence line from record: " + name);

        line = reader.readLine();
        if(line==null)
        	throw new RuntimeException("Missing comment line from record: " + name);

        if (line.charAt(0)=='+') {
            comment = line.substring(1);
        } else {
            throw new RuntimeException("Invalid FASTQ comment line: " + line);
        }

        quality = reader.readLine();
        if(quality==null)
        	throw new RuntimeException("Missing quality line from record: " + name);

        current = new FastqRecord(name, sequence, comment, quality, phredOffset);
    }

    public int getProgress() {
    	if(atEOF.get())
    		return 100;
    	
    	long bytesRead=posTrackInputStream.getPosition();
    	
    	return (int)(((float) bytesRead / fileLength) * 100);    
    }

    
    private void accumulateHistogram(FastqRecord rec)
    {
    	int quals[]=rec.getQualityAsInteger(false);
    	
    	for(int i: quals)
    		qualHistogram[i]++;
    }
    
    public int determinePhredOffset()
    {
    	int phred33Total=0;
    	int phred64Total=0;

    	for(int i=33;i<=58;i++)
    		phred33Total+=qualHistogram[i];
    	
    	for(int i=80;i<=104;i++)
    		phred64Total+=qualHistogram[i];
    	
    	if(phred33Total==0 && phred64Total>0)
    		return 64;

    	if(phred64Total==0 && phred33Total>0)
    		return 33;
    	
    	return 0;
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
        
        if(phredOffset==0)
        	{
        	deque.clear();
        	qualHistogram=new int[256];
        	
        	for(int i=0;i<PREREAD_COUNT;i++)
        		{
        		parseOne();
        		if(current!=null)
        			{
        			deque.add(current);
        			accumulateHistogram(current);
        			}
        		}
        	}
        parseOne();
    }

    public void close() throws IOException {
        reader.close();
    }

    public boolean hasNext() {
        return (!deque.isEmpty()) || (current != null);
    }

    public FastqRecord next() throws IOException {
    	if(deque.isEmpty())
    		{
    		FastqRecord current = this.current;
    		parseOne();

    		return current;
    		}
    	else
    		{
    		FastqRecord rec=deque.poll();
    		
    		if(rec!=null)
    			rec.setPhredOffset(phredOffset);
    		
    		return rec;
    		}
    }

}
