package org.umcn.me.output.vcf;

import java.util.Locale;

public class MobsterRecordVCFWrapper implements VCFFunctions {

	private MobsterRecord record;
	private String identifier = ".";
	
	public final static String VCFHEADER = "##fileformat=VCFv4.1\n" + 
			"##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n" + 
			"##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n" + 
			"##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n" + 
			"##ALT=<ID=INS:ME:HERV,Description=\"Insertion of HERV element\">\n" + 
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n" + 
			"##INFO=<ID=SAMPS,Number=1,Type=String,Description=\"Number of supporting reads per sample. Sample names and read counts are seperated by a colon (:) while samples are seperated by |\">\n" + 
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n" + 
			"##INFO=<ID=SUP,Number=4,Type=Integer,Description=\"Number of supporting reads in the form of SUPPORTING READS ON 5' SIDE (BOTH DISCORDANT AND SPLIT), SUPPORTING READS ON 3' SIDE (BOTH DISCORDANT AND SPLIT), SUPPORTING SPLIT READS ON 5' SIDE, SUPPORTING SPLIT READS ON 3' SIDE\">\n" + 
			"##INFO=<ID=POLYA,Number=1,Type=Integer,Description=\"Number of split reads containing a POLYA stretch\">\n" + 
			"##INFO=<ID=TSD,Number=1,Type=String,Description=\"Whether it is suspected that a target site duplication or deletion is present\">\n" + 
			"##INFO=<ID=CLLEN,Number=2,Type=Integer,Description=\"Length of the 5' side cluster composed of discordant anchors, length of the 3' side cluster composed of discordant anchors. Value is 0 when there are no discordant anchors on the particular side of the insertion\">\n" + 
			"##INFO=<ID=ORIGIN,Number=3,Type=Integer,Description=\"Specification of how the discordant pairs supporting the event were aligned in the original BAM file: Number of pairs both aligning uniquely but discordantly, number of pairs where one end aligns multiple times to the reference, number of pairs where one end is unmapped\">\n" + 
			"##INFO=<ID=CLIPPED,Number=6,Type=Float,Description=\"Information about the clipped reads in the form of: The maximum distance between clipping positions of reads which are clipped on left side, The maximum distance between clipping positions of reads clipped on the right side, Fraction of left-clipped reads with same clipping position, Fraction of right-clipped reads with same clipping position, Mean base quality of clipped subsequences, Mean length of clipped part of read\">\n" +
			"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
	
	
	public MobsterRecordVCFWrapper(MobsterRecord record){
		if (record.getCorrectedBorder5() > record.getBorder3()){
			System.err.println("[WARN] Corrected 5' CIPOS position is bigger than 3' CIPOS position for: " + record);
		}
		this.record = record;
	
	}

	@Override
	public String getChromosome() {
		return record.getChromosome();
	}

	@Override
	public Integer getPosition() {
		return record.getInsertionPoint();
	}

	@Override
	public String getID() {
		return this.identifier;
	}

	@Override
	public String getReference() {
		return ".";
	}

	@Override
	public String getAlternative() {
		return "<INS:ME:" + record.getMobileElement() + ">";
	}

	@Override
	public String getQuality() {
		return ".";
	}

	@Override
	public String getFilter() {
		return "PASS";
	}

	@Override
	public String getInfo() {
		return VCFDefinitions.IMPRECISE + VCFDefinitions.INFO_SEPERATOR + this.getSampleCounts() + VCFDefinitions.INFO_SEPERATOR + 
				this.getCIPos() + VCFDefinitions.INFO_SEPERATOR + 
				this.getSupportingReads() + VCFDefinitions.INFO_SEPERATOR + this.getPolyASupport() +
				VCFDefinitions.INFO_SEPERATOR + this.getTSD() + VCFDefinitions.INFO_SEPERATOR + this.getClusterLengths() +
				VCFDefinitions.INFO_SEPERATOR + this.getOrigins() + VCFDefinitions.INFO_SEPERATOR + this.getClippingInfo();
	}
	
	
	public String getFirstSevenFields(){
	   return this.getChromosome() + "\t" + this.getPosition() + "\t" + this.getID() + "\t" + this.getReference() + 
			   "\t" + this.getAlternative() + "\t" + this.getQuality() + "\t" + this.getFilter();
   }
	   
		   
	//Check below for INFO functions
	  
	public String getCIPos(){		
		return VCFDefinitions.CIPOS + "=" + Integer.toString(this.record.getCorrectedBorder5()) + "," + Integer.toString(this.record.getBorder3());
	}
	
	public String getSupportingReads(){
		return VCFDefinitions.SUPPORTING_READS + "=" + Integer.toString(this.record.getCluster5Hits()) + "," + Integer.toString(this.record.getCluster3Hits())
				+ "," + Integer.toString(this.record.getSplit5Hits()) + "," + Integer.toString(this.record.getSplit3Hits());
	}
	
	public String getPolyASupport(){
		return VCFDefinitions.POLY_A + "=" + Integer.toString(this.record.getTotalPolyAReads());
	}
	
	public String getOrigins(){
		return VCFDefinitions.READ_ORIGIN + "=" + Integer.toString(this.record.getuUpairs()) + "," + Integer.toString(this.record.getuMpairs())
				+ "," + Integer.toString(this.record.getuXpairs());
	}
	
	public String getSampleCounts(){
		return VCFDefinitions.SAMPLECOUNTS + "=" + this.record.getCorrectedSampleCounts();
	}
	
	public String getClippingInfo(){
		return VCFDefinitions.CLIPPED_INFO + "=" + Integer.toString(this.record.getLeftClippedMaxDist()) + "," +
				Integer.toString(this.record.getRightClippedMaxDist()) + "," + this.record.getLeftClippedSamePos()
				+ "," + this.record.getRightClippedSamePos() + "," + String.format( Locale.ENGLISH, "%.2f", this.record.getClippedAvgQual())
				+ "," + String.format( Locale.ENGLISH, "%.2f", this.record.getClippedAvgLength());
	}
	
	/**
	 * First get cluster5 length then get cluster3 length
	 * @return
	 */
	public String getClusterLengths(){
		return VCFDefinitions.CLUSTER_LENGTH + "=" + Integer.toString(this.record.getCluster5Length()) + "," + Integer.toString(this.record.getCluster3Length());
	}
	
	public String getTSD(){
		return VCFDefinitions.TSD + "=" + this.record.getTsd();
	}

   
   public String toString(){
	   return this.getFirstSevenFields() + "\t" + this.getInfo();
   }

	
}
