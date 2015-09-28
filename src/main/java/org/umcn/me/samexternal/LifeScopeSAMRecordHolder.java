package org.umcn.me.samexternal;



import net.sf.samtools.SAMRecord;

/**
 * @author Djie
 */

public class LifeScopeSAMRecordHolder extends NrMappingsSAMRecordHolder {

	private int nr_mappings;
	private String nr_mappings_attribute = SAMDefinitions.LIFESCOPE_NRHITS_ATTRIBUTE;
	
	public LifeScopeSAMRecordHolder(SAMRecord record){
		super(record);
		parseNumberOfMappings();
	}
	
	private void parseNumberOfMappings(){
		
		if (this.sam_record.getReadUnmappedFlag()){
			this.nr_mappings = 0;
		}else{
			this.nr_mappings = Integer.parseInt(this.sam_record.getAttribute(nr_mappings_attribute).toString());
		}
	}
	
	@Override
	public boolean isMappedUniquely() {
		return (nr_mappings == 1);
	}

	@Override
	public boolean isMappedMultiple() {
		return (nr_mappings > 1);
	}

	@Override
	public int getNumberOfMappings() {
		return nr_mappings;
	}

}
