package org.umcn.me.output.vcf;

import org.beanio.annotation.Field;
import org.beanio.annotation.Record;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Container for a Mobster line
 * @author Djie Tjwan Thung
 *
 */

@Record(minOccurs=1)
public class MobsterRecord implements Comparable<MobsterRecord> {
	
	
//0		1				2					3		4		5		6		7				8				9				10		
//Chr	Mobile Element	Insert Point	border5	border3	merged	sample	sample_counts	cluster5 length	cluster3 length	cluster5 hits
//11			12				13			14			15			16			17			18							19					20									
//cluster3 hits	split5 hits	split3 hits	polyA5 hits	polyT5 hits	polyA3 hits	polyT3 hits	original discordant unique	original multiple	original unmapped
//21						22					23		
//leftclipped max dist	rightclipped max dist	leftclipped same pos
//24						25						26				27		
//rightclipped same pos	clipped avg qual	clipped avg length	target site duplication
	
	private String chromosome;	
	private String mobileElement;	
	private int insertionPoint;	
	private int border5;	
	private int border3;
	private int endInsertionPoint;
	private int endBorder5;
	private int endBorder3;
	private String sample;
	private String sampleCounts;
	private String nonSupportingSampleCounts;
	private String cluster5Length;	
	private String cluster3Length;	
	private String cluster5Hits;	
	private String cluster3Hits;	
	private int split5Hits;
	private int split3Hits;
	private int polyA5hits;
	private int polyT5hits;
	private int polyA3hits;
	private int polyT3hits;
	private int uUpairs;
	private int uMpairs;
	private int uXpairs;

	private int leftClippedMaxDist;
	private int rightClippedMaxDist;
	private String leftClippedSamePos;
	private String rightClippedSamePos;
	
	private double clippedAvgQual;
	private double clippedAvgLength;
	private String tsd;
	private int tsdLen;
	private String tsdSeq;
	private String sampleVAFs;
	private String somatic;

	public static MobsterRecord fromVCF(String line, String samples) {
		MobsterRecord record = new MobsterRecord();

		record.setChromosome(getTag("CHROM", line));

		int insPoint = Integer.parseInt(getTag("POS", line));
		int endInsPoint = Integer.parseInt(getTag("END", line));
		record.setInsertionPoint(insPoint);
		record.setEndInsertionPoint(endInsPoint);

		record.setBorder5(insPoint + Integer.parseInt(getTag("CIPOS", line).split(",")[0]));
		record.setBorder3(insPoint + Integer.parseInt(getTag("CIPOS", line).split(",")[1]));
		record.setEndBorder5(insPoint + Integer.parseInt(getTag("CIEND", line).split(",")[0]));
		record.setEndBorder3(insPoint + Integer.parseInt(getTag("CIEND", line).split(",")[1]));

		record.setMobileElement(getTag("ID", line).replace("<INS:ME:", "").replace(">",""));

		record.setCluster5Length(getTag("CLLEN", line).split(",")[0]);
		record.setCluster3Length(getTag("CLLEN", line).split(",")[1]);

		record.setCluster5Hits(getTag("SUP", line).split(",")[0]);
		record.setCluster3Hits(getTag("SUP", line).split(",")[1]);
		record.setSplit5Hits(Integer.parseInt(getTag("SUP", line).split(",")[2]));
		record.setSplit3Hits(Integer.parseInt(getTag("SUP", line).split(",")[3]));

		record.setPolyA5hits(Integer.parseInt(getTag("POLYA", line).split(",")[0]));
		record.setPolyT5hits(Integer.parseInt(getTag("POLYA", line).split(",")[1]));
		record.setPolyA3hits(Integer.parseInt(getTag("SUP", line).split(",")[2]));
		record.setPolyT5hits(Integer.parseInt(getTag("SUP", line).split(",")[3]));

		record.setTSD(getTag("TSD", line));
		record.setTSDlen(Integer.parseInt(getTag("TSDLEN", line)));
		record.setTSDseq(getTag("TSDSEQ", line));

		record.setuUpairs(Integer.parseInt(getTag("ORIGIN", line).split(",")[0]));
		record.setuXpairs(Integer.parseInt(getTag("ORIGIN", line).split(",")[1]));
		record.setuMpairs(Integer.parseInt(getTag("ORIGIN", line).split(",")[2]));

		record.setLeftClippedMaxDist(Integer.parseInt(getTag("CLIPPED", line).split(",")[0]));
		record.setRightClippedMaxDist(Integer.parseInt(getTag("CLIPPED", line).split(",")[1]));
		record.setLeftClippedSamePos(getTag("CLIPPED", line).split(",")[2]);
		record.setRightClippedSamePos(getTag("CLIPPED", line).split(",")[3]);
		record.setClippedAvgQual(Double.parseDouble(getTag("CLIPPED", line).split(",")[4]));
		record.setClippedAvgLength(Double.parseDouble(getTag("CLIPPED", line).split(",")[5]));

		record.setSample(samples);
		Pattern pattern = Pattern.compile("./.:[^\t\n]+");
		Matcher matcher = pattern.matcher(line);
		ArrayList<String> sampleCounts = new ArrayList<String>();
		ArrayList<String> sampleVAFs = new ArrayList<String>();
		while(matcher.find()){
			String genotype = matcher.group(0);
			sampleCounts.add(genotype.split(":")[1]);
			sampleVAFs.add(genotype.split(":")[2]);
		}
		record.setSampleCountsFromVCF(sampleCounts.toArray(new String[0]));
		record.setSampleVAFsFromVCF(sampleVAFs.toArray(new String[0]));

		return record;
	}

	public static String getTag(String tag, String line){
		String[] splitLine = line.split("\t");
		switch(tag){
			case "CHROM":
				return splitLine[0];
			case "POS":
				return splitLine[1];
			case "ID":
				return splitLine[2];
			case "REF":
				return splitLine[3];
			case "ALT":
				return splitLine[4];
			case "QUAL":
				return splitLine[5];
			case "FILTER":
				return splitLine[6];
			default:
				Pattern pattern = Pattern.compile("(^|;)"+tag+"=([^;]+)");
				Matcher matcher = pattern.matcher(splitLine[7]);
				matcher.find();
				return matcher.group(2).trim();
		}

	}

	public String getChromosome() { return chromosome; }
	
	@Field(at=0, required=true)
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	
	public String getMobileElement() { return mobileElement; }
	
	@Field(at=1, required=true)
	public void setMobileElement(String mobileElement) {
		this.mobileElement = mobileElement;
	}
	
	/**
	 * Get 1-based insertion coordinate
	 * if 1-based insertion coordinate is < 1, correct to 1
	 * @return
	 */
	public int getInsertionPoint() {
		if (this.insertionPoint < 1){
			return 1;
		}
		return this.insertionPoint;
	}
	
	@Field(at=2, required=true)
	public void setInsertionPoint(int insertionPoint) {
		this.insertionPoint = insertionPoint;
	}
	
	public int getBorder5() {
		return border5;
	}
	
	@Field(at=3, required=true)
	public void setBorder5(int border5) {
		this.border5 = border5;
	}
	
	public int getBorder3() {
		return border3;
	}

	@Field(at=4, required=true)
	public void setBorder3(int border3) {
		this.border3 = border3;
	}

	public int getEndInsertionPoint() {
		if (this.endInsertionPoint < 1){
			return 1;
		}
		return this.endInsertionPoint;
	}

	@Field(at=5, required=true)
	public void setEndInsertionPoint(int endInsertionPoint) {
		this.endInsertionPoint = endInsertionPoint;
	}

	public int getEndBorder5() {
		return endBorder5;
	}

	@Field(at=6, required=true)
	public void setEndBorder5(int endBorder5) {
		this.endBorder5 = endBorder5;
	}

	public int getEndBorder3() {
		return endBorder3;
	}

	@Field(at=7, required=true)
	public void setEndBorder3(int endBorder3) {
		this.endBorder3 = endBorder3;
	}

	public String getSample() { return sample; }
	
	@Field(at=9, required=true)
	public void setSample(String sample) {
		this.sample = sample;
	}

	public String getSampleCounts() { return sampleCounts; }

	@Field(at=10, required=true)
	public void setSampleCounts(String sampleCounts) {
		ArrayList<String> sampleCountsArray = new ArrayList<>();
		for(String sample: sampleCounts.split(", "))
			sampleCountsArray.add(sample.split("=")[1]);
		this.sampleCounts = String.join(", ", sampleCountsArray);
	}
	public void setSampleCountsFromVCF(String[] sampleCounts) {
		this.sampleCounts = String.join(", ", sampleCounts);
	}

	public String getNonSupportingSampleCounts() { return nonSupportingSampleCounts; }

	@Field(at=11, required=true)
	public void setNonSupportingSampleCounts(String nonSupportingSampleCounts) {
		ArrayList<String> nonSupportingSampleCountsArray = new ArrayList<>();
		for(String sample: nonSupportingSampleCounts.split(", "))
			nonSupportingSampleCountsArray.add(sample.split("=")[1]);
		this.nonSupportingSampleCounts = String.join(", ", nonSupportingSampleCountsArray);
	}
	public void setNonSupportingSampleCountsFromVCF(String[] nonSupportingSampleCounts) {
		this.nonSupportingSampleCounts = String.join(", ", nonSupportingSampleCounts);
	}


	public int getCluster5Length() {
		return NAConverter.toInt(this.cluster5Length);
	}
	
	@Field(at=12, required=true)
	public void setCluster5Length(String cluster5Length) {
		this.cluster5Length = cluster5Length;
	}
	
	public int getCluster3Length() {
		return NAConverter.toInt(this.cluster3Length);
	}
	
	@Field(at=13, required=true)
	public void setCluster3Length(String cluster3Length) {
		this.cluster3Length = cluster3Length;
	}
	
	public int getCluster5Hits() {
		return NAConverter.toInt(this.cluster5Hits);
	}
	
	@Field(at=14, required=true)
	public void setCluster5Hits(String cluster5Hits) {
		this.cluster5Hits = cluster5Hits;
	}
	
	public int getCluster3Hits() {
		return NAConverter.toInt(this.cluster3Hits);
	}
	
	@Field(at=15, required=true)
	public void setCluster3Hits(String cluster3Hits) {
		this.cluster3Hits = cluster3Hits;
	}
	
	public int getSplit5Hits() {
		return split5Hits;
	}
	
	@Field(at=16, required=true)
	public void setSplit5Hits(int split5Hits) {
		this.split5Hits = split5Hits;
	}
	
	public int getSplit3Hits() {
		return split3Hits;
	}
	
	@Field(at=17, required=true)
	public void setSplit3Hits(int split3Hits) {
		this.split3Hits = split3Hits;
	}
	
	public int getPolyA5hits() {
		return polyA5hits;
	}
	
	@Field(at=18, required=true)
	public void setPolyA5hits(int polyA5hits) {
		this.polyA5hits = polyA5hits;
	}
	
	public int getPolyT5hits() {
		return polyT5hits;
	}
	
	@Field(at=19, required=true)
	public void setPolyT5hits(int polyT5hits) {
		this.polyT5hits = polyT5hits;
	}
	
	public int getPolyA3hits() {
		return polyA3hits;
	}
	
	@Field(at=20, required=true)
	public void setPolyA3hits(int polyA3hits) {
		this.polyA3hits = polyA3hits;
	}
	
	public int getPolyT3hits() {
		return polyT3hits;
	}
	
	@Field(at=21, required=true)
	public void setPolyT3hits(int polyT3hits) {
		this.polyT3hits = polyT3hits;
	}
	
	public int getuUpairs() {
		return uUpairs;
	}
	
	@Field(at=22, required=true)
	public void setuUpairs(int uUpairs) {
		this.uUpairs = uUpairs;
	}
	
	public int getuMpairs() {
		return uMpairs;
	}
	
	@Field(at=23, required=true)
	public void setuMpairs(int uMpairs) {
		this.uMpairs = uMpairs;
	}
	
	public int getuXpairs() {
		return uXpairs;
	}
	
	@Field(at=24, required=true)
	public void setuXpairs(int uXpairs) {
		this.uXpairs = uXpairs;
	}
	
	public int getLeftClippedMaxDist() {
		return leftClippedMaxDist;
	}

	@Field(at=25, required=true)
	public void setLeftClippedMaxDist(int leftClippedMaxDist) {
		this.leftClippedMaxDist = leftClippedMaxDist;
	}

	public int getRightClippedMaxDist() {
		return rightClippedMaxDist;
	}

	@Field(at=26, required=true)
	public void setRightClippedMaxDist(int rightClippedMaxDist) {
		this.rightClippedMaxDist = rightClippedMaxDist;
	}

	public String getLeftClippedSamePos() {
		return leftClippedSamePos;
	}

	@Field(at=27, required=true)
	public void setLeftClippedSamePos(String leftClippedSamePos) {
		this.leftClippedSamePos = leftClippedSamePos;
	}

	public String getRightClippedSamePos() {
		return rightClippedSamePos;
	}

	@Field(at=28, required=true)
	public void setRightClippedSamePos(String rightClippedSamePos) {
		this.rightClippedSamePos = rightClippedSamePos;
	}
	
	
	public double getClippedAvgQual() {
		return clippedAvgQual;
	}
	
	@Field(at=29, required=true)
	public void setClippedAvgQual(double clippedAvgQual) {
		this.clippedAvgQual = clippedAvgQual;
	}
	public double getClippedAvgLength() {
		return clippedAvgLength;
	}
	
	@Field(at=30, required=true)
	public void setClippedAvgLength(double clippedAvgLength) {
		this.clippedAvgLength = clippedAvgLength;
	}
	
	public String getTSD() { return tsd; }

	@Field(at=31, required=true)
	public void setTSD(String tsd) { this.tsd = tsd; }

	public int getTSDlen() { return tsdLen; }

	@Field(at=32, required=true)
	public void setTSDlen(int tsdLen) { this.tsdLen = tsdLen; }

	public String getTSDseq() { return tsdSeq; }

	@Field(at=33, required=true)
	public void setTSDseq(String tsdSeq) { this.tsdSeq = tsdSeq; }

	public String getSampleVAFs() { return sampleVAFs; }

	@Field(at=34, required=true)
	public void setSampleVAFs(String sampleVAFs) {
		ArrayList<String> sampleVAFsArray = new ArrayList<>();
		for(String sample: sampleVAFs.split(", "))
			sampleVAFsArray.add(sample.split("=")[1]);
		this.sampleVAFs = String.join(", ", sampleVAFsArray);
	}
	public void setSampleVAFsFromVCF(String[] vafs) { this.sampleVAFs = String.join(", ", vafs); }

	public String getSomatic() { return somatic; }

	@Field(at=35, required=false, minOccurs = 0)
	public void setSomatic(String somatic) { this.somatic = somatic; }

	//Other methods
	
	/**
	 * Correct negative integer border5's
	 * @return
	 */
	public int getCorrectedBorder5(){
		int newBorder5 = this.getBorder5();
		if (newBorder5 <= 0){
			newBorder5 = 1;
		}
		return newBorder5;
	}
	public int getCorrectedEndBorder5(){
		int newEndBorder5 = this.getEndBorder5();
		if (newEndBorder5 <= 0){
			newEndBorder5 = 1;
		}
		return newEndBorder5;
	}
	
	/**
	 * Sample counts with legal values for VCF format
	 * @return
	 */
//	public String getCorrectedSampleCounts(){
//		return this.sampleCounts.replaceAll("=", ":").replaceAll(", ", "|");
//	}
	
	public int getTotalPolyAReads(){
		return this.getPolyA3hits() + this.getPolyA5hits() + this.getPolyT3hits() + this.getPolyT5hits();
	}

	
	public String toString(){
		return this.chromosome + "\t" + this.border5 + "\t" +  this.border3 + "\t" +  this.mobileElement;
	}

	@Override
	public int compareTo(MobsterRecord rec2) {
		//returns negative if "this" object is less than "that" object
	    //returns 0 if they are equal
	    //returns positive if "this" object is greater than "that" object
		
		if (rec2.getChromosome().equals(this.getChromosome())){
			return this.getInsertionPoint() - rec2.getInsertionPoint();
		}else{
			return this.getChromosome().compareTo(rec2.getChromosome());
		}
	}
}
