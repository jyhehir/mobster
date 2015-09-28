package org.umcn.me.util;

import java.security.InvalidParameterException;
import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;

import org.umcn.me.samexternal.SAMSilentReader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.picard.sam.SamFileHeaderMerger;


public class BAMCollection {

	List<BAMSample> bams = new ArrayList<BAMSample>();
	
	
	public BAMCollection(String[] files, String[] names) throws InvalidParameterException, FileNotFoundException {
		
		this(files, names, false);
		
	}
	
	public BAMCollection(String[] files, String[] names, boolean prefix) throws InvalidParameterException, FileNotFoundException {
				
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
			BAMSample sample = new BAMSample(file, name);
			
			if (prefix) sample.setPrefixReadGroupId(Integer.toString(i) + ".");
			bams.add(sample);
		}
	}
	
	/**
	 * 
	 * @param order
	 * @return SAMFileHeader
	 */
	public SAMFileHeader getMergedHeader(SAMFileHeader.SortOrder order) {
		
		List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();		
		
		for (BAMSample sample : this.bams){
			List<SAMReadGroupRecord> newReadGroups = new ArrayList<SAMReadGroupRecord>();
			SAMSilentReader reader = new SAMSilentReader(sample.bam);
			SAMFileHeader header = reader.getFileHeader();
			
			//If no readgroup is present, make an artificial one
			if (header.getReadGroups().isEmpty()){
				SAMReadGroupRecord record = new SAMReadGroupRecord (sample.getPrefixReadGroupId());
				record.setSample(sample.sample);
				newReadGroups.add(record);
			}else{
				for (SAMReadGroupRecord rg : header.getReadGroups()){
					rg.setSample(sample.sample);
					SAMReadGroupRecord rgDup = new SAMReadGroupRecord (sample.getPrefixReadGroupId() + rg.getReadGroupId(), rg);
					newReadGroups.add(rgDup);
				}
			}
			reader.close();		
			
			header.setReadGroups(newReadGroups);
			headers.add(header);

		}
		
		SamFileHeaderMerger merger = new SamFileHeaderMerger(order, headers, true);//true: do merge the sequence dictionaries
		return merger.getMergedHeader();
	}
	
	public String getPrefixReadGroupIdFromBam(BAMSample bam){
		
		int index = this.bams.indexOf(bam);
		
		return this.bams.get(index).getPrefixReadGroupId();
		
	}
	
	public List<BAMSample> getCloneOfBAMSampleList(){
		
		List<BAMSample> copyList = new ArrayList<BAMSample>();
		
		for (BAMSample sample : this.bams){
			copyList.add(new BAMSample(sample));
		}
		
		return copyList;
		
	}
	
	
}
