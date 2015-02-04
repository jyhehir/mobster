package org.umcn.me.util;

import java.security.InvalidParameterException;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;

import org.umcn.gen.sam.SAMSilentReader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.picard.sam.SamFileHeaderMerger;


public class BAMCollection {

	List<BAMSample> bams = new ArrayList<BAMSample>();
	
	public BAMCollection(String[] files, String[] names) throws InvalidParameterException, FileNotFoundException {
	
		if (files.length != names.length){
			throw new InvalidParameterException("Number of files does not equal number of sample names");
		}
		
		for (int i=0; i < files.length; i++){
			
			String fileName = files[i];
			File file = new File(fileName);
			String name = names[i];
			
			if (fileName.isEmpty() || name.isEmpty()){
				throw new InvalidParameterException("Name of file or sample may not be empty");
			}
			if (!file.exists()){
				throw new FileNotFoundException("File does not exist: " + file);
			}
			
			bams.add(new BAMSample(file, name));
		}
		
	}
	
	public SAMFileHeader getMergedHeader(SAMFileHeader.SortOrder order){
		
		List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();		
		
		for (BAMSample sample : this.bams){
			SAMSilentReader reader = new SAMSilentReader(sample.bam);
			SAMFileHeader header = reader.getFileHeader();
			
			for (SAMReadGroupRecord rg : header.getReadGroups()){
				rg.setSample(sample.sample);
			}
			headers.add(header);
			
			reader.close();
		}
		
		SamFileHeaderMerger merger = new SamFileHeaderMerger(order, headers, true);//true: do merge the sequence dictionaries
		return merger.getMergedHeader();
	}
	
}
