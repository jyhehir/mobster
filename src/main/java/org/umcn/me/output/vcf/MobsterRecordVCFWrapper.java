package org.umcn.me.output.vcf;

import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.util.InvalidNucleotideSequenceException;
import org.umcn.me.util.ReferenceGenome;

import java.util.Arrays;
import java.util.Locale;

public class MobsterRecordVCFWrapper implements VCFFunctions {

	private MobsterRecord record;
	private String identifier = ".";
	private ReferenceGenome referenceGenome;
	private String[] allSamples;

	public MobsterRecordVCFWrapper(MobsterRecord record){
		if (record.getCorrectedBorder5() > record.getBorder3()){
			System.err.println("[WARN] Corrected 5' CIPOS position is bigger than 3' CIPOS position for: " + record);
		}
		if (record.getCorrectedEndBorder5() > record.getEndBorder3()){
			System.err.println("[WARN] Corrected 5' CIEND position is bigger than 3' CIPEND position for: " + record);
		}
		this.record = record;
	}

	public MobsterRecordVCFWrapper(MobsterRecord record, String[] allSamples){
		this(record);
		this.allSamples = allSamples;
	}

	public MobsterRecordVCFWrapper(MobsterRecord record, String[] allSamples, ReferenceGenome referenceGenome){
		this(record, allSamples);
		this.referenceGenome = referenceGenome;
	}

	public static String getHeader(){
		return VCFDefinitions.VCFHEADER;
	}
	public static String getHeader(String[] samples){
		return VCFDefinitions.VCFHEADER + String.join("\t", samples) + "\n";
	}

	public String getFirstSevenFields(){
		return this.getChromosome() + "\t" + this.getPosition() + "\t" + this.getID() + "\t" + this.getReference() +
				"\t" + this.getAlternative() + "\t" + this.getQuality() + "\t" + this.getFilter();
	}

	@Override
	public String getInfo() {
		return (isSomatic()? (VCFDefinitions.SOMATIC + VCFDefinitions.INFO_SEPERATOR) : "") + (isImprecise()? (VCFDefinitions.IMPRECISE + VCFDefinitions.INFO_SEPERATOR) : "") +
				this.getSVtype() + VCFDefinitions.INFO_SEPERATOR + this.getCIPOS() + VCFDefinitions.INFO_SEPERATOR +
				this.getEnd() + VCFDefinitions.INFO_SEPERATOR + this.getCIEND() + VCFDefinitions.INFO_SEPERATOR +
				this.getSupportingReads() + VCFDefinitions.INFO_SEPERATOR + this.getPolyASupport() + VCFDefinitions.INFO_SEPERATOR +
				this.getTSD() + VCFDefinitions.INFO_SEPERATOR + this.getTSDlen() + VCFDefinitions.INFO_SEPERATOR +
				this.getTSDseq() + VCFDefinitions.INFO_SEPERATOR + this.getClusterLengths() + VCFDefinitions.INFO_SEPERATOR +
				this.getOrigins() + VCFDefinitions.INFO_SEPERATOR + this.getClippingInfo();
	}

	public String getGenotype(){
		StringBuilder genotypeColumn = new StringBuilder("GT:AD:AF");
		String[] sampleNames = this.record.getSample().split(", ");
		String[] sampleCounts = this.record.getSampleCounts().split(", ");
		String[] nonSupportingSampleCounts = this.record.getNonSupportingSampleCounts().split(", ");
		String[] sampleVAFs = this.record.getSampleVAFs().split(", ");

		for(String sampleName: this.allSamples){
			int sampleIndex = Arrays.asList(sampleNames).indexOf(sampleName);
			String sampleGenotype = "./.";
			String sampleCount = ".";
			String nonSupportingSampleCount = ".";
			String sampleVAF = ".";
			if(sampleIndex != -1){
				sampleGenotype = "0/1";
				sampleCount = sampleCounts[sampleIndex];
				nonSupportingSampleCount = nonSupportingSampleCounts[sampleIndex];
				sampleVAF  = sampleVAFs[sampleIndex];
			}

			genotypeColumn.append("\t").append(sampleGenotype).append(":").append(nonSupportingSampleCount).append(",").append(sampleCount).append(":");
			if(!sampleVAF.equals("-1"))
				genotypeColumn.append(sampleVAF);
			else
				genotypeColumn.append(".");
		}

		return genotypeColumn.toString();
	}

	public String getSample(){
		return this.record.getSample();
	}

	//Methods for the primary fields
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
		if(this.referenceGenome == null){
			return ".";
		}
		else {
			try {
				return Character.toString(this.referenceGenome.getBaseAt(getChromosome(), getPosition()));
			} catch (Exception e) {
				System.out.println(e.toString());
				return ".";
			}
		}
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
		if(this.record.getSomatic() == null)
			return "PASS";
		else
			switch(this.record.getSomatic()){
				case MobilePrediction.GERMLINE:
					return "germline";
				case MobilePrediction.ARTIFACT:
					return "normal_artifact";
				case MobilePrediction.SOMATIC:
					return "PASS";
				default:
					return "PASS";
		}
	}
		   
	//Check below for INFO functions
	public boolean isSomatic(){
		return MobilePrediction.SOMATIC.equals(this.record.getSomatic()) || MobilePrediction.ARTIFACT.equals(this.record.getSomatic());
	}

	public boolean isImprecise(){
		return this.record.getInsertionPoint() != this.record.getCorrectedBorder5() || this.record.getInsertionPoint() != this.record.getBorder3()
			|| this.record.getEndInsertionPoint() != this.record.getCorrectedEndBorder5() || this.record.getEndInsertionPoint() != this.record.getEndBorder3();
	}

	public String getSVtype(){
		return VCFDefinitions.SVTYPE + "=INS:ME:" + record.getMobileElement();
	}

	public String getCIPOS(){
		return VCFDefinitions.CIPOS + "=" + Integer.toString(this.record.getCorrectedBorder5() - this.record.getInsertionPoint()) + "," + Integer.toString(this.record.getBorder3()  - this.record.getInsertionPoint());
	}

	public String getEnd(){
		return VCFDefinitions.END + "=" + Integer.toString(this.record.getEndInsertionPoint());
	}

	public String getCIEND(){
		return VCFDefinitions.CIEND + "=" + Integer.toString(this.record.getCorrectedEndBorder5() - this.record.getEndInsertionPoint()) + "," + Integer.toString(this.record.getEndBorder3() - this.record.getEndInsertionPoint());
	}
	
	public String getSupportingReads(){
		return VCFDefinitions.SUPPORTING_READS + "=" + Integer.toString(this.record.getCluster5Hits()) + "," + Integer.toString(this.record.getCluster3Hits())
				+ "," + Integer.toString(this.record.getSplit5Hits()) + "," + Integer.toString(this.record.getSplit3Hits());
	}
	
	public String getPolyASupport(){
		return VCFDefinitions.POLY_A + "=" + Integer.toString(this.record.getPolyA5hits()) + "," + Integer.toString(this.record.getPolyT5hits())
				+ "," + Integer.toString(this.record.getPolyA3hits()) + "," + Integer.toString(this.record.getPolyT3hits());
	}
	
	public String getOrigins(){
		return VCFDefinitions.READ_ORIGIN + "=" + Integer.toString(this.record.getuUpairs()) + "," + Integer.toString(this.record.getuMpairs())
				+ "," + Integer.toString(this.record.getuXpairs());
	}

//	public String getSampleCounts(){
//		return VCFDefinitions.SAMPLECOUNTS + "=" + this.record.getCorrectedSampleCounts();
//	}

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
		return VCFDefinitions.TSD + "=" + this.record.getTSD();
	}

	public String getTSDlen(){ return VCFDefinitions.TSDLEN + "=" + this.record.getTSDlen(); }

	public String getTSDseq(){ return VCFDefinitions.TSDSEQ + "=" + this.record.getTSDseq(); }

//	public String getSampleVAFs(){
//		return VCFDefinitions.VAF + "=" + this.record.getSampleVAFs();
//	}

   public String toString(){
	   return this.getFirstSevenFields() + "\t" + this.getInfo() + "\t" + this.getGenotype();
   }
}
