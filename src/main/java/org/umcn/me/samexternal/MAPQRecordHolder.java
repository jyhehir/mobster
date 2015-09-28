package org.umcn.me.samexternal;

import net.sf.samtools.SAMRecord;


/**
 * Define uniqueness / multi mapping on mapq
 * @author Djie Tjwan Thung
 *
 */
public class MAPQRecordHolder extends NrMappingsSAMRecordHolder {

	private int nr_mappings;
	
	/**
	 * MAPQ threshold: >= threshold everything is considered  unique
	 * Below this threshold everything is considered non-unique
	 * 
	 */
	private int mapq_threshold = 20;
	
	public MAPQRecordHolder(SAMRecord samRecord){
		super(samRecord);
		parseNumberOfMappings();
	}
	
	public MAPQRecordHolder(SAMRecord samRecord, int minClipping, int maxClipping){
		super(samRecord, minClipping, maxClipping);
		parseNumberOfMappings();
	}
	
	/**
	 * 
	 * @param samRecord
	 * @param minClipping
	 * @param maxClipping
	 * @param mapqUnique everything >= threshold is considered unique; everything < is considered multiple
	 */
	public MAPQRecordHolder(SAMRecord samRecord, int minClipping, int maxClipping, int mapqUnique){
		super(samRecord, minClipping, maxClipping);
		this.mapq_threshold = mapqUnique;
		parseNumberOfMappings();
		
	}
	
	private void parseNumberOfMappings(){

		if (this.sam_record.getReadUnmappedFlag()){
			this.nr_mappings = 0;
		}else{
			if (this.sam_record.getMappingQuality() >= this.mapq_threshold){
				this.nr_mappings = 1;
			}else {
				this.nr_mappings = 999;
			}
		}
	}
	
	@Override
	public boolean isMappedUniquely() {
		if (this.nr_mappings == 1){
			return true;
		}
		return false;
	}
	

	@Override
	public boolean isMappedMultiple() {
		if (this.nr_mappings > 1){
			return true;
		}
		return false;
	}

	@Override
	public int getNumberOfMappings() {
		return this.nr_mappings;
	}
}
