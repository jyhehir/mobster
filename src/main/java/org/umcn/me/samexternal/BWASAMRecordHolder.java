package org.umcn.me.samexternal;


import net.sf.samtools.SAMRecord;
/**
 * @author Djie
 */

public class BWASAMRecordHolder extends NrMappingsSAMRecordHolder {
	
	private String mapping_info_attribute = SAMDefinitions.BWA_NR_OPTIMAL_HITS_ATTRIBUTE;
	private int nr_mappings;
	
	public BWASAMRecordHolder(SAMRecord samRecord){
		super(samRecord);
		parseNumberOfMappings();
	}
	
	public BWASAMRecordHolder(SAMRecord samRecord, int minClipping, int maxClipping){
		super(samRecord, minClipping, maxClipping);
		parseNumberOfMappings();
	}
	
	private void parseNumberOfMappings(){

		if (this.sam_record.getReadUnmappedFlag()){
			this.nr_mappings = 0;
		}else{
			if(this.sam_record.getAttribute(this.mapping_info_attribute) != null){
				this.nr_mappings = Integer.parseInt(
						this.sam_record.getAttribute(this.mapping_info_attribute).toString());
			}else{
				this.nr_mappings = 1;
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
