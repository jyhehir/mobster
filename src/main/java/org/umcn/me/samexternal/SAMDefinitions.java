package org.umcn.me.samexternal;

/**
 * @author Djie
 */
public class SAMDefinitions {

	public static final String MAPPING_TOOL_BWA = "bwa";
	public static final String MAPPING_TOOL_MOSAIK = "mosaik";
	public static final String MAPPING_TOOL_LIFESCOPE = "lifescope";
	public static final String MAPPING_TOOL_UNSPECIFIED = "unspecified";
	
	public static final String[] MAPPING_TOOLS = {MAPPING_TOOL_MOSAIK, MAPPING_TOOL_BWA, MAPPING_TOOL_LIFESCOPE, MAPPING_TOOL_UNSPECIFIED};
	public static final String MOSAIK_MAPPINGINFO_ATTRIBUTE = "ZA";
	public static final String LIFESCOPE_NRHITS_ATTRIBUTE = "NH";
	public static final String BWA_NR_OPTIMAL_HITS_ATTRIBUTE = "X0";
	public static final String BWA_NR_SUBOPTIMAL_HITS_ATTRIBUTE = "X1";
	public static final int MOSAIK_MAPPINGINFO_READ1 = 0;
	public static final int MOSAIK_MAPPINGINFO_READ2 = 1;
	public static final int MOSAIK_MAPPINGINFO_NRHITS = 4;
	public static final String READ_NUMBER_SEPERATOR = "-";
	public static final String UNIQUE_UNIQUE_MAPPING = "UU";
	public static final String UNIQUE_UNMAPPED_MAPPING = "UX";
	public static final String UNIQUE_MULTIPLE_MAPPING = "UM";
	public static final String SPLIT_MAPPING = "S";
	public static final char LEFT_CLIPPED = 'L';
	public static final char RIGHT_CLIPPED = 'R';

}
