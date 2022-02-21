package org.umcn.me.util;


public abstract class MobileDefinitions {
	
	public static final String[] MOBILE_CATEGORIES = {"ALU", "SVA", "L1", "HERV"};
	
	public static final String OTHER_MOBILE_CATEGORY_ATTRIBUTE = ".ME_CAT.";
	
	public static final String SAM_TAG_MOBILE_HIT = "MH";
	public static final String SAM_TAG_CLUSTER_LENGTH = "CL";
	public static final String SAM_TAG_CLUSTER_HITS = "CH";
	public static final String SAM_TAG_MOBILE = "ME";
	public static final String SAM_TAG_SPLIT_CLUSTER = "YC";
	public static final String SAM_TAG_UNIQUE_HITS = "YN";
	public static final String SAM_TAG_MULTIPLE_HITS = "YO";
	public static final String SAM_TAG_UNMAPPED_HITS = "YP";
	
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_HITS =  "YA";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_HITS = "YB";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_MAX_DISTANCE = "YD";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_MAX_DISTANCE = "YF";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_FRAC_SAME_DISTANCE = "YG";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_FRAC_SAME_DISTANCE = "YH";
	public static final String SAM_TAG_SPLIT_CLIPPED_AVG_QUAL = "YI";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_ENDS = "YJ";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_ENDS = "YK";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_MEDIAN_END = "YL";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_MEDIAN_END = "YM";
	public static final String SAM_TAG_SPLIT_AVG_CLIPPED_LEN = "YQ";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_POLYT = "YR";
	public static final String SAM_TAG_SPLIT_LEFTCLIPPED_POLYA = "YS";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_POLYT = "YT";
	public static final String SAM_TAG_SPLIT_RIGHTCLIPPED_POLYA = "YU";
	public static final String SAM_TAG_SAMPLENAME = "YY";
	public static final String SAM_TAG_SAMPLECOUNT = "YX";
	public static final String SAM_TAG_READNAMES = "YZ";
	
	//properties Mobster
	public static final String READ_LENGTH = "READ_LENGTH";
	public static final String USE_READ_LENGTH = "USE_READ_LENGTH";
	
	//properties PotentialMEIReadFinder
	
	public static final String INFILE = "IN_FILE";
	public static final String RESOURCES_DIR = "RESOURCES_DIR";
	public static final String OUTFILE = "OUT_FILE";
	public static final String VCF = "VCF_OUT";
	public static final String MAPPING_TOOL = "MAPPING_TOOL";
	public static final String USE_SPLIT = "USE_SPLIT";
	public static final String TMP = "TEMP_DIR";
	public static final String MIN_CLIPPING = "MINIMUM_CLIP_LENGTH";
	public static final String MAX_CLIPPING = "MAXIMUM_CLIP_LENGTH";
	public static final String MEMORY = "MAX_RECORDS_IN_RAM";
	public static final String MIN_QUAL = "MINIMUM_AVG_QUALITY";
	public static final String MIN_MAPQ_ANCHOR = "MINIMUM_MAPQ_ANCHOR";
	
	public static final String QUERY_SORT_INPUT = "QUERY_SORT_INPUT";
	
	
	//properties RefAndMEPairFinder
	
//	File single = null;
//	File multiple = null;
//	File filtered = null;
//	String out;
//	String sampleName = "";
//	String tool = "";
//	Boolean paired;
	
	public static final String INFILE_FROM_MOBIOME_MAPPING = "BAM_FROM_MOBIOME_MAPPING";
	public static final String INFILE_FROM_POTENTIAL_MEI_FINDER = "BAM_FROM_POTENTIALMEIFINDER";
	public static final String SAMPLE_NAME = "SAMPLENAME";
	public static final String MOBIOME_MAPPING_TOOL = "MOBIOME_MAPPING_TOOL";
	public static final String PAIRED_END = "PAIRED_END";
	public static final String POLY_A_LENGTH = "MINIMUM_POLYA_LENGTH";
	public static final String POLY_A_MAX_MISMATCHES = "MAXIMUM_MISMATCHES_POLYA";
	
	//properties AnchorClusterer
	//clusterBam = OUT_FILE
	//sample already in sample_name property
	public static final String ANCHOR_FILE = "ANCHOR_BAM_FILE";
	public static final String SPLIT_ANCHOR_FILE = "ANCHOR_SPLIT_BAM_FILE";
	public static final String READS_PER_CLUSTER = "READS_PER_CLUSTER";
	public static final String DISCORDANT_OVERLAP = "DISCORDANT_CLUSTERS_MAX_OVERLAP";
	public static final String DISCORDANT_DISTANCE = "DISCORDANT_CLUSTER_MAX_DISTANCE";
	public static final String MEAN_FRAGMENT_LENGTH = "MEAN_FRAGMENT_LENGTH";
	public static final String SD_FRAGMENT_LENGTH = "SD_FRAGMENT_LENGTH";
	public static final String LENGTH_99PROCENT_OF_FRAGMENTS = "LENGTH_99PROCENT_OF_FRAGMENTS";
	public static final String MAX_SPACING_CLIPPED_READS = "MAX_SPACING_OF_CLIPPED_READS";
	public static final String MAX_OVERLAP_CLIPPED_CLUSTERS = "MAX_OVERLAP_OF_CLIPPED_CLUSTERS";
	public static final String MAX_DISTANCE_CLIPPED_CLUSTERS = "MAX_DISTANCE_OF_CLIPPED_CLUSTERS";
	public static final String SEARCH_AREA = "NEIGHBORHOOD_WINDOW_BP";
	public static final String MINIMUM_TOTAL_HITS = "MINIMUM_SUPPORTING_READS";
	public static final String MINIMUM_INITIAL_SPLIT_CLUSTER_READS = "MINIMUM_INITIAL_SPLIT_CLUSTER_READS";
	public static final String REPEAT_MASK_FILE = "REPEATMASK_FILE";
        public static final String PREFIX_REFERENCE = "PREFIX_REFERENCE";
	
	public static final String MOBIOME_MAPPING_CMD = "MOBIOME_MAPPING_CMD";
	public static final String PICARD_COLLECT_INSERT_METRICS = "PICARD_COLLECT_INSERT_SIZE_METRICS_JAR";
	public static final String USE_PICARD = "USE_PICARD";

	public static final String MULTIPLE_SAMPLE_CALLING = "MULTIPLE_SAMPLE_CALLING";

	public static final String MULTIPLE_SAMPLE_CALLING_STRINGENT = "MULTIPLE_SAMPLE_CALLING_STRINGENT";
	
	public static final String DEFAULT_SEP=",";
	
	
	//properties for GRIPS
	
	public static final String GRIPS_MIN_PERCENT_SAME_REF_MATE_MAPPING = "GRIPS_MIN_PERC_SAME_REF_MAPPING";
	public static final String GRIPS_DISCARD_UNIQUE_MULTIPLE = "GRIPS_NO_UM_PAIRS";
	public static final String GRIPS_MAX_REFSEQ_MAPPING = "GRIPS_MAX_REFSEQ";
	public static final String GRIPS_MAX_SOURCE_GENES = "GRIPS_MAX_SOURCE_GENES";
	public static final String GRIPS_MAX_DIFF_MATE_MAPPINGS = "GRIPS_MAX_DIFF_MATE_MAPPINGS";
	public static final String FILTER_OVERLAPPING_PREDICTIONS = "FILTER_OVERLAPPING_PREDICTIONS";
	public static final String GRIPS_MODE = "GRIPS_MODE";
	
	public static final String GRIPS_SELFCHAIN = "SELFCHAIN";
	public static final String GRIPS_TRANSCRIPT = "TRANSCRIPT";
	public static final String GRIPS_BLACKLIST = "BLACKLIST";
	public static final String GRIPS_REPMASK = "GRIPS_REPMASK";
	public static final String GRIPS_SKIP_UU_EXCLUSION = "GRIPS_SKIP_UU_EXCLUSION"; 
	
	public static final String CLEANUP_FILES = "CLEANUP";

}
