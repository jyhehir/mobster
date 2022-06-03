package org.umcn.me.output.vcf;

import org.beanio.annotation.Field;
import org.beanio.annotation.Record;

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
	private String sample;	
	private String sampleCounts;	
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
	private double vaf;

	public String getChromosome() { return chromosome; }
	
	@Field(at=0, required=true)
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	
	public String getMobileElement() {
		return mobileElement;
	}
	
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
		this.insertionPoint = insertionPoint + 1;
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
	
	public String getSample() {
		return sample;
	}
	
	@Field(at=6, required=true)
	public void setSample(String sample) {
		this.sample = sample;
	}
	
	public String getSampleCounts() {
		return sampleCounts;
	}
	
	@Field(at=7, required=true)
	public void setSampleCounts(String sampleCounts) {
		this.sampleCounts = sampleCounts;
	}
	
	public int getCluster5Length() {
		return NAConverter.toInt(this.cluster5Length);
	}
	
	@Field(at=8, required=true)
	public void setCluster5Length(String cluster5Length) {
		this.cluster5Length = cluster5Length;
	}
	
	public int getCluster3Length() {
		return NAConverter.toInt(this.cluster3Length);
	}
	
	@Field(at=9, required=true)
	public void setCluster3Length(String cluster3Length) {
		this.cluster3Length = cluster3Length;
	}
	
	public int getCluster5Hits() {
		return NAConverter.toInt(this.cluster5Hits);
	}
	
	@Field(at=10, required=true)
	public void setCluster5Hits(String cluster5Hits) {
		this.cluster5Hits = cluster5Hits;
	}
	
	public int getCluster3Hits() {
		return NAConverter.toInt(this.cluster3Hits);
	}
	
	@Field(at=11, required=true)
	public void setCluster3Hits(String cluster3Hits) {
		this.cluster3Hits = cluster3Hits;
	}
	
	public int getSplit5Hits() {
		return split5Hits;
	}
	
	@Field(at=12, required=true)
	public void setSplit5Hits(int split5Hits) {
		this.split5Hits = split5Hits;
	}
	
	public int getSplit3Hits() {
		return split3Hits;
	}
	
	@Field(at=13, required=true)
	public void setSplit3Hits(int split3Hits) {
		this.split3Hits = split3Hits;
	}
	
	public int getPolyA5hits() {
		return polyA5hits;
	}
	
	@Field(at=14, required=true)
	public void setPolyA5hits(int polyA5hits) {
		this.polyA5hits = polyA5hits;
	}
	
	public int getPolyT5hits() {
		return polyT5hits;
	}
	
	@Field(at=15, required=true)
	public void setPolyT5hits(int polyT5hits) {
		this.polyT5hits = polyT5hits;
	}
	
	public int getPolyA3hits() {
		return polyA3hits;
	}
	
	@Field(at=16, required=true)
	public void setPolyA3hits(int polyA3hits) {
		this.polyA3hits = polyA3hits;
	}
	
	public int getPolyT3hits() {
		return polyT3hits;
	}
	
	@Field(at=17, required=true)
	public void setPolyT3hits(int polyT3hits) {
		this.polyT3hits = polyT3hits;
	}
	
	public int getuUpairs() {
		return uUpairs;
	}
	
	@Field(at=18, required=true)
	public void setuUpairs(int uUpairs) {
		this.uUpairs = uUpairs;
	}
	
	public int getuMpairs() {
		return uMpairs;
	}
	
	@Field(at=19, required=true)
	public void setuMpairs(int uMpairs) {
		this.uMpairs = uMpairs;
	}
	
	public int getuXpairs() {
		return uXpairs;
	}
	
	@Field(at=20, required=true)
	public void setuXpairs(int uXpairs) {
		this.uXpairs = uXpairs;
	}
	
	public int getLeftClippedMaxDist() {
		return leftClippedMaxDist;
	}

	@Field(at=21, required=true)
	public void setLeftClippedMaxDist(int leftClippedMaxDist) {
		this.leftClippedMaxDist = leftClippedMaxDist;
	}

	public int getRightClippedMaxDist() {
		return rightClippedMaxDist;
	}

	@Field(at=22, required=true)
	public void setRightClippedMaxDist(int rightClippedMaxDist) {
		this.rightClippedMaxDist = rightClippedMaxDist;
	}

	public String getLeftClippedSamePos() {
		return leftClippedSamePos;
	}

	@Field(at=23, required=true)
	public void setLeftClippedSamePos(String leftClippedSamePos) {
		this.leftClippedSamePos = leftClippedSamePos;
	}

	public String getRightClippedSamePos() {
		return rightClippedSamePos;
	}

	@Field(at=24, required=true)
	public void setRightClippedSamePos(String rightClippedSamePos) {
		this.rightClippedSamePos = rightClippedSamePos;
	}
	
	
	public double getClippedAvgQual() {
		return clippedAvgQual;
	}
	
	@Field(at=25, required=true)
	public void setClippedAvgQual(double clippedAvgQual) {
		this.clippedAvgQual = clippedAvgQual;
	}
	public double getClippedAvgLength() {
		return clippedAvgLength;
	}
	
	@Field(at=26, required=true)
	public void setClippedAvgLength(double clippedAvgLength) {
		this.clippedAvgLength = clippedAvgLength;
	}
	
	public String getTsd() { return tsd; }

	@Field(at=27, required=true)
	public void setTsd(String tsd) { this.tsd = tsd; }

	public double getVAF() { return vaf; }

	@Field(at=28, required=true)
	public void setVAF(double vaf) { this.vaf = vaf; }

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
	
	/**
	 * Sample counts with legal values for VCF format
	 * @return
	 */
	public String getCorrectedSampleCounts(){
		return this.sampleCounts.replaceAll("=", ":").replaceAll(", ", "|");
	}
	
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
