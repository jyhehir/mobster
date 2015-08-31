package org.umcn.me.samexternal;

import java.io.File;

import net.sf.samtools.SAMFileReader;

public class SAMSilentReader extends SAMFileReader {
	
	public SAMSilentReader(File samOrBam){
		super(samOrBam);
		this.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	}
	
	public SAMSilentReader(File bam, File index){
		super(bam, index);
		this.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	}
	
}
