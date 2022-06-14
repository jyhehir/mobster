package org.umcn.me.output.vcf;

public class VCFDefinitions {

	public static final String INFO_SEPERATOR = ";";
	
	/**
	 * Header of first 9 fields
	 */
	public static final String BASIC_HEADER = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	
	public static final String CIPOS = "CIPOS";
	public static final String END = "END";
	public static final String CIEND = "CIEND";
	public static final String IMPRECISE = "IMPRECISE";
	
	//MEI INFO fields
	public static final String SVTYPE = "SVTYPE";
	public static final String CLUSTER_LENGTH = "CLLEN";
	public static final String SUPPORTING_READS = "SUP";
	public static final String POLY_A="POLYA";
	public static final String READ_ORIGIN="ORIGIN";
	public static final String CLIPPED_INFO="CLIPPED";
	public static final String TSD = "TSD";
	public static final String TSDLEN = "TSDLEN";
	public static final String TSDSEQ = "TSDSEQ";
	public static final String SAMPLECOUNTS = "SAMPS";
	public static final String VAF = "AF";


}
