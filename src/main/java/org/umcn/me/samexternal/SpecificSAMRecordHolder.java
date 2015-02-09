package org.umcn.me.samexternal;

import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import net.sf.samtools.SAMRecord;

/**
 * @author Djie
 */

public class SpecificSAMRecordHolder {
	protected SAMRecord sam_record;
	protected int soft_clipping_left;
	protected int soft_clipping_right;
	
	public SpecificSAMRecordHolder(SAMRecord samRecord){
		this.sam_record = samRecord;
		this.soft_clipping_left = getNumberSoftClippedBases(true);
		this.soft_clipping_right = getNumberSoftClippedBases(false);
	}
	
	public SAMRecord getSAMRecord(){
		return this.sam_record;
	}
	
	public int getSoftClippingLeft(){
		return this.soft_clipping_left;
	}
	
	public int getSoftClippingRight(){
		return this.soft_clipping_right;		
	}
	
	public boolean isLargetsClippedSeqLeft(){
		return (this.soft_clipping_left >= this.soft_clipping_right);
	}
	
	public String getSoftClippedSeqLeft(){
		return this.sam_record.getReadString().substring(0, this.soft_clipping_left);
	}
	
	public String getSoftClippedSeqRight(){
		String readString = this.sam_record.getReadString();
		int readLength = readString.length();
		return readString.substring(readLength - this.soft_clipping_right);
	}
	
	public String getSoftClippedQualLeft(){
		return this.sam_record.getBaseQualityString().substring(0, this.soft_clipping_left);
	}
	
	public String getSoftClippedQualRight(){
		String qualString = this.sam_record.getBaseQualityString();
		int qualLength = qualString.length();
		return qualString.substring(qualLength - this.soft_clipping_right);
	}
	
	public String getLargestClippedSeq(){
		if (this.soft_clipping_left >= this.soft_clipping_right){
			return this.getSoftClippedSeqLeft();
		}else if(this.soft_clipping_right > this.soft_clipping_left){
			return this.getSoftClippedSeqRight();
		}
		return "";
	}
	
	public String getLargestClippedQual(){
		if (this.soft_clipping_left >= this.soft_clipping_right){
			return this.getSoftClippedQualLeft();
		}else if(this.soft_clipping_right > this.soft_clipping_left){
			return this.getSoftClippedQualRight();
		}
		return "";
	}
	
	private int getNumberSoftClippedBases(boolean left){
		Pattern softClippedPattern;
		Matcher match;
		String matchedString;
		int softClippedNr = 0;
		
		String cigar = this.sam_record.getCigarString();
		
		if(left){
			softClippedPattern = Pattern.compile("^\\d+S");
		}else{
			softClippedPattern = Pattern.compile("\\d+S$");
		}
		
		match = softClippedPattern.matcher(cigar);
		
		if (match.find()){
			matchedString = match.group();
			softClippedNr = Integer.parseInt(matchedString.substring(0, 
					matchedString.length() - 1)); //remove trailing S
		}
		
		return softClippedNr;
			
	}
	
	@Override
	public String toString(){
		return this.sam_record.getSAMString();
	}
	
	/**
	 * New functions for NIPD below:
	 */
	
	public boolean isOnReference(String ref){
		return this.sam_record.getReferenceName().toString().equals(ref) && this.isMapped();
	}
	
	public boolean isInReferenceList(Vector<String> references){
		
		String recordReference = this.sam_record.getReferenceName();
		if(references.contains(recordReference)){
			return true;
		}else{
			for(String region : references){
				if(region.contains(":")){
					String[] firstSplit = region.split(":");
					String ref = firstSplit[0];
					int start = Integer.parseInt(firstSplit[1].split("-")[0]);
					int end = Integer.parseInt(firstSplit[1].split("-")[1]);
					
					if(isOnReferenceAndWithinBounds(ref, start, end)){
						return true;
					}
				}
			}
		}
		return false;
		
	}
	
	public boolean isOnReferenceAndWithinBounds(String ref, int start, int end){
		if(isOnReference(ref)){
			int recordStart = this.sam_record.getAlignmentStart();
			int recordEnd = this.sam_record.getAlignmentEnd();
			
			//
			//q  ----		  ----- (both queries will return true
			//ref	-----------
			if (recordEnd >= start && recordStart <= end){
				return true;
			}
		}
		return false;
	}
	
	public boolean isMapped(){
		return !this.sam_record.getReadUnmappedFlag();
	}
	
	
	/**
	 * Note when record has no NM attribute this function will return true!
	 * @param maxEdit
	 * @return
	 */
	public boolean hasLowerEditDistanceThan(int maxEdit){
		if(this.sam_record.getAttribute("NM") == null){
			return true;
		}else{
			return (Integer.parseInt(this.sam_record.getAttribute("NM").toString()) <= maxEdit);
		}
	}
	
	public boolean hasXbasesHigherThanQualY(int xBases, int yQual){
		Vector<Integer> quals = QualityProcessing.decodeQualStringToVector(this.sam_record.getBaseQualityString(), 33);
		return (QualityProcessing.passesQualityCheck(quals, xBases, yQual));
	}
	
	public boolean hasMAPQHigherThan(int mapq){
		return (this.sam_record.getMappingQuality() >= mapq);
	}
	
	public boolean hasMeanQualHigherThan(double threshold){
		Vector<Integer> quals = QualityProcessing.decodeQualStringToVector(this.sam_record.getBaseQualityString(), 33);
		return (QualityProcessing.calculateAverage(quals) >= threshold);
	}
	
}
