package org.umcn.me.output.vcf;

public class VCFDefinitions {

	public final static String VCFHEADER = "##fileformat=VCFv4.2\n" +
			"##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element.\">\n" +
			"##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element.\">\n" +
			"##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element.\">\n" +
			"##ALT=<ID=INS:ME:HERV,Description=\"Insertion of HERV element.\">\n" +
			"##FILTER=<ID=PASS,Description=\"When combined tumor/normal calling is enabled, only insertions with supporting reads in the tumor and none in the normal sample will get a PASS.\">\n" +
			"##FILTER=<ID=germline,Description=\"The number of reads supporting the insertion in the normal sample is higher than the minimum and is therefore likely present in the germline.\">\n" +
			"##FILTER=<ID=normal_artifact,Description=\"The number of reads supporting the insertion in the normal sample is higher than 0 but lower than the minimum required to set it as germline (non-somatic). A possible explanation for these reads is that tumor contamination occurred in the normal sample or another artifact.\">\n" +
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation.\">\n" +
			"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"A somatic variation.\">\n" +
			"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" +
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants with no supporting split reads for this position.\">\n" +
			"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant. If it greater than POS, a target site duplication or deletion is present for the sequence from POS+1 until this coordinate.\">\n" +
			"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants with no supporting split reads for this position.\">\n" +
			"##INFO=<ID=SUP,Number=4,Type=Integer,Description=\"Number of supporting reads in the form of: SUPPORTING READS ON 5' SIDE (BOTH DISCORDANT AND SPLIT), SUPPORTING READS ON 3' SIDE (BOTH DISCORDANT AND SPLIT), SUPPORTING SPLIT READS ON 5' SIDE, SUPPORTING SPLIT READS ON 3' SIDE\">\n" +
			"##INFO=<ID=POLYA,Number=4,Type=Integer,Description=\"Number of split reads containing a POLYA/T stretch in the form of: POLYA ON 5' SIDE, POLYT ON 5' SIDE, POLYA ON 3' SIDE, POLYT ON 3' SIDE\">\n" +
			"##INFO=<ID=TSD,Number=1,Type=String,Description=\"Whether it is suspected that a target site duplication or deletion is present flanking the insertion.\">\n" +
			"##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description=\"The length of the target site duplication or deletion if present.\">\n" +
			"##INFO=<ID=TSDSEQ,Number=1,Type=String,Description=\"The sequence of the target site duplication or deletion if present.\">\n" +
			"##INFO=<ID=CLLEN,Number=2,Type=Integer,Description=\"Length of the 5' side cluster composed of discordant anchors, length of the 3' side cluster composed of discordant anchors. Value is 0 when there are no discordant anchors on the particular side of the insertion.\">\n" +
			"##INFO=<ID=ORIGIN,Number=3,Type=Integer,Description=\"Specification of how the discordant pairs supporting the event were aligned in the original BAM file: Number of pairs both aligning uniquely but discordantly, number of pairs where one end aligns multiple times to the reference, number of pairs where one end is unmapped.\">\n" +
			"##INFO=<ID=CLIPPED,Number=6,Type=Float,Description=\"Information about the clipped reads in the form of: The maximum distance between clipping positions of reads which are clipped on left side, The maximum distance between clipping positions of reads clipped on the right side, Fraction of left-clipped reads with same clipping position, Fraction of right-clipped reads with same clipping position, Mean base quality of clipped subsequences, Mean length of clipped part of read.\">\n" +
			"##FORMAT=<ID=GT,Type=String,Number=1,Type=String,Description=\"Genotype, not called\">\n" +
			"##FORMAT=<ID=AD,Type=Integer, Number=1,Type=String,Description=\"The number of reads spanning the clusters that were not supporting the insertion, followed by the number of supporting reads per sample\">\n" +
			"##FORMAT=<ID=AF,Type=Float, Number=1,Type=Float,Description=\"The allele frequency based on raw read counts per sample. Likely an underestimation.\">\n" +
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\t";

	public static final String INFO_SEPERATOR = ";";

	/**
	 * Header of first 9 fields
	 */
	public static final String BASIC_HEADER = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	
	public static final String CIPOS = "CIPOS";
	public static final String END = "END";
	public static final String CIEND = "CIEND";
	public static final String IMPRECISE = "IMPRECISE";
	public static final String SOMATIC = "SOMATIC";

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
