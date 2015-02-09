package org.umcn.me.samexternal;


/**
 * @author Djie
 */

import net.sf.samtools.SAMRecord;

public class MosaikSAMRecordHolder extends NrMappingsSAMRecordHolder {
	
	private String mapping_info_attribute = SAMDefinitions.MOSAIK_MAPPINGINFO_ATTRIBUTE;
	private String mapping_info_attribute_value = null;;
	private int nr_mappings;
	
	public MosaikSAMRecordHolder(SAMRecord samRecord){
		super(samRecord);
		this.setMappingInfoAttributeValue();
		this.parseNrMappingsTag();
	}
	
	public MosaikSAMRecordHolder(SAMRecord samRecord, int minClipping, int maxClipping){
		super(samRecord, minClipping, maxClipping);
		this.setMappingInfoAttributeValue();
		this.parseNrMappingsTag();
	}
	
	private void parseNrMappingsTag(){
		String readNrMappingInfo;
		if (this.mapping_info_attribute_value != null) {
			if (this.sam_record.getReadPairedFlag() && this.sam_record.getFirstOfPairFlag()) {
				readNrMappingInfo = this.mapping_info_attribute_value.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ1];
			} else if (this.sam_record.getReadPairedFlag() && this.sam_record.getSecondOfPairFlag()) {
				readNrMappingInfo = this.mapping_info_attribute_value.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ2];
			} else{ //If single-end fragment
				readNrMappingInfo = this.mapping_info_attribute_value;
			}
			this.nr_mappings = Integer
					.parseInt(readNrMappingInfo.split(";")[SAMDefinitions.MOSAIK_MAPPINGINFO_NRHITS]);
		}else{ //No Mapping Info Attribute on Mosaik SAM Record means no mapping
			this.nr_mappings = 0;
		}	
	}
	
	private void setMappingInfoAttributeValue(){
		Object mappingInfoAttribute = this.sam_record.getAttribute(this.mapping_info_attribute);
		if (mappingInfoAttribute != null){
			this.mapping_info_attribute_value = mappingInfoAttribute.toString();
		}
	}
	
	@Override
	public boolean isMappedUniquely() {
		if(this.nr_mappings == 1 && !this.sam_record.getReadUnmappedFlag()){
			return true;
		}
		return false;
	}

	@Override
	public boolean isMappedMultiple() {
		if(this.nr_mappings > 1 && !this.sam_record.getReadUnmappedFlag()){
			return true;
		}
		return false;
	}

	@Override
	public int getNumberOfMappings() {
		return this.nr_mappings;
	}

}
