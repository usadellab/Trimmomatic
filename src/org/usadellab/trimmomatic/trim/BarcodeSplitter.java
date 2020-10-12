/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.usadellab.trimmomatic.trim;

import java.util.HashMap;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 *
 * @author marc
 */
public class BarcodeSplitter extends AbstractSingleRecordTrimmer {

    private HashMap<String, String> barcodes;
    private int maxMisMatch = 0;
    private boolean clipOffBarcodes = true;    

    public BarcodeSplitter(HashMap<String, String> barcodes, int mism, boolean clip) {
        this.barcodes = barcodes;
        this.maxMisMatch = mism;
        this.clipOffBarcodes = clip;
    }
    
    public HashMap<String, String> getBarcodeMap() {
        return barcodes;
    }

    @Override
    public FastqRecord processRecord(FastqRecord entry) {
        
        String seq = entry.getSequence();
        
        for (String label : barcodes.keySet()) {
            int codeLength = barcodes.get(label).length();
            if (maxMisMatch == 0) {
                if (entry.getSequence().startsWith(barcodes.get(label))) {
                    
                    if (clipOffBarcodes) {
                    FastqRecord outentry = new FastqRecord(entry, codeLength, seq.length() - codeLength);
                        outentry.setBarcodeLabel(label);                    
                        return outentry;
                    } else {
                        entry.setBarcodeLabel(label);                    
                        return entry;
                    }
                }
            } else {
                if (getMisMatches(barcodes.get(label), entry.getSequence(), barcodes.get(label).length()) <= maxMisMatch) {
                    
                    if (clipOffBarcodes) {
                    FastqRecord outentry = new FastqRecord(entry, codeLength, seq.length() - codeLength);
                        outentry.setBarcodeLabel(label);                    
                        return outentry;
                    } else {
                        entry.setBarcodeLabel(label);                    
                        return entry;
                    }
                }
            }
        }        
        // entry did not match any barcodes - write to UNKNOWN
        entry.setBarcodeLabel("UNKNOWN");
        return entry;
    }

    private int getMisMatches(String a, String b, int range) {
        int mismatch = 0;
        for (int i = 0; i < range; i++) {
            mismatch += (a.charAt(i) == b.charAt(i)) ? 0 : 1;
        }
        return mismatch;
    }
}
