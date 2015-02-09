package org.umcn.me.samexternal;

import net.sf.samtools.SAMRecord;

/**
 * @author Djie
 */

public class NrMappingsSAMRecordHolderFactory {
	
	public static NrMappingsSAMRecordHolder makeNrMappingsSAMRecordHolder(SAMRecord record,
			String mappingTool) throws UnknownParamException{
		
		if (SAMDefinitions.MAPPING_TOOL_MOSAIK.equals(mappingTool)){
			return new MosaikSAMRecordHolder(record);
		} else if(SAMDefinitions.MAPPING_TOOL_BWA.equals(mappingTool)){
			return new BWASAMRecordHolder(record);
		} else if (SAMDefinitions.MAPPING_TOOL_LIFESCOPE.equals(mappingTool)){
			return new LifeScopeSAMRecordHolder(record);
		} else if (SAMDefinitions.MAPPING_TOOL_UNSPECIFIED.equals(mappingTool)){
			return new MAPQRecordHolder(record);
		}
		//add below support for new subclasses of NrMappingsSAMRecordHolder
		
		throw new UnknownParamException("No support for the following mapping tool: " + mappingTool);
	}
	
	public static NrMappingsSAMRecordHolder makeNrMappingsSAMRecordHolder(SAMRecord record,
			String mappingTool, int minClipping, int maxClipping) throws UnknownParamException{
		
		if (SAMDefinitions.MAPPING_TOOL_MOSAIK.equals(mappingTool)){
			return new MosaikSAMRecordHolder(record, minClipping, maxClipping);
		} else if(SAMDefinitions.MAPPING_TOOL_BWA.equals(mappingTool)){
			return new BWASAMRecordHolder(record, minClipping, maxClipping);
		} else if (SAMDefinitions.MAPPING_TOOL_UNSPECIFIED.equals(mappingTool)){
			return new MAPQRecordHolder(record, minClipping, maxClipping);
		}
		
		//add below support for new subclasses of NrMappingsSAMRecordHolder
		
		throw new UnknownParamException("No support for the following mapping tool: " + mappingTool);
	}
	
	//TODO: Dirty method needs fixing: Method is implemented so a factory method can 
	//be made initializing a MAPQRecordHolder.
	//This is dirty as the mapq parameter is not passed on to MOSAIK or bwa... i.e. it doesnt do anything,
	//while the user may expect this.
	public static NrMappingsSAMRecordHolder makeNrMappingsSAMRecordHolder(SAMRecord record,
			String mappingTool, int minClipping, int maxClipping, int mapq) throws UnknownParamException{
		
		if (SAMDefinitions.MAPPING_TOOL_UNSPECIFIED.equals(mappingTool)){
			return new MAPQRecordHolder(record, minClipping, maxClipping, mapq);
		}else{
			return makeNrMappingsSAMRecordHolder(record, mappingTool, minClipping, maxClipping);
		}
		
	}
	
	
}
