package org.umcn.me.splitread;

import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.umcn.me.sam.MobileSAMTag;
import org.umcn.me.samexternal.QualityProcessing;
import org.umcn.me.util.MobileDefinitions;

import net.sf.samtools.SAMRecord;

/**
 * This class represents all reads which are partially hard-clipped in bam files.
 * 
 * 
 * @author Djie Thung
 *
 */
public class ClippedRead  {

	private boolean start_clipped;
	private int clipped_length;
	private String quality_string_encoded;
	private Vector<Integer> quality_vector_decoded;
	private int phred_offset = 33;
	private SAMRecord rec = null;
	private boolean color_space;
	
	public ClippedRead(SAMRecord rec, boolean softClipped, boolean colorSpace) throws InvalidHardClipCigarException {
		this.rec = rec;
		this.color_space = colorSpace;
		processCigarAndQual(softClipped);
	}
	
	public ClippedRead(SAMRecord rec, boolean softClipped, boolean colorSpace, int offset) throws InvalidHardClipCigarException {
		this(rec, softClipped, colorSpace);
		this.phred_offset = offset;
	}
	
	private void processCigarAndQual(Boolean softClipped) throws InvalidHardClipCigarException{
		processClippedSideAndLength(softClipped);
		processQual();
	}
	
	private int processClippedLength(boolean softClipped, boolean left){
		Pattern clippedPattern;
		Matcher match;
		String matchedString;
		StringBuilder pattern = new StringBuilder();
		int clippedNr = 0;
		
		String cigar = this.rec.getCigarString();
		
		if(left){
			pattern.append("^");
		}
		
		if (softClipped){
			pattern.append("\\d+S");
		}else{
			pattern.append("\\d+H");
		}

		if(!left){
			pattern.append("$");
		}

		
		clippedPattern = Pattern.compile(pattern.toString());
		match = clippedPattern.matcher(cigar);
		
		if (match.find()){
			matchedString = match.group();
			clippedNr = Integer.parseInt(matchedString.substring(0, 
					matchedString.length() - 1)); //remove trailing S
		}
		
		return clippedNr;
	}
	
	
	/**
	 * Determine which side of the read is clipped and how
	 * long the clipping sequence is.
	 * If both sides are clipped used the longest clipped end.
	 * 
	 * The obtained information is used to set the variables
	 * start_clipped and clipped_length
	 * @throws InvalidHardClipCigarException
	 */
	private void processClippedSideAndLength(boolean softClipped) throws InvalidHardClipCigarException{
		
		int clippedLengthStart = 0;
		int clippedLengthEnd = 0;
		
		clippedLengthStart = processClippedLength(softClipped, true);
		clippedLengthEnd = processClippedLength(softClipped, false);
		
		if(clippedLengthStart == 0 && clippedLengthEnd == 0){
			throw new InvalidHardClipCigarException("Cigar should contain at least one hardclip segment.");
		}

		if (clippedLengthStart >= clippedLengthEnd){
			this.clipped_length = clippedLengthStart;
			this.start_clipped = true;
		}else{
			this.clipped_length = clippedLengthEnd;
			this.start_clipped = false;
		}						
	}
	
	/**
	 * Set the quality attributes for this class to both the encoded string
	 * and a vector of PHRED integer quality scores based on the phred_offset
	 * default is 33 (ASCII value of 33 == PHRED score of 0), but this value
	 * can be changed through the constructor
	 */
	private void processQual(){
		if(this.color_space){
			processQualColorSpace();
		}else{
			processQualBaseSpace();
		}
	}

	private void processQualBaseSpace(){
		if(this.start_clipped){
			this.quality_string_encoded = this.rec.getBaseQualityString().substring(0, this.clipped_length);
			this.quality_vector_decoded = QualityProcessing.decodeQualStringToVector(this.quality_string_encoded, phred_offset);
		}else{
			String qualString = this.rec.getBaseQualityString();
			int qualLength = qualString.length();
			this.quality_string_encoded = qualString.substring(qualLength - this.clipped_length);
			this.quality_vector_decoded = QualityProcessing.decodeQualStringToVector(this.quality_string_encoded, phred_offset);
		}
	}
	
	private void processQualColorSpace() {
		boolean negStrand = this.rec.getReadNegativeStrandFlag();
		
		if ((negStrand && !this.start_clipped) || (!negStrand && this.start_clipped)){
			//Do not use first base as it contains the contaminated base next to the primer
			this.quality_string_encoded = this.rec.getAttribute("CQ").toString().substring(1, this.clipped_length);
			this.quality_vector_decoded = QualityProcessing.decodeQualStringToVector(this.quality_string_encoded, phred_offset);			
		}else if((negStrand && this.start_clipped) || (!negStrand && !this.start_clipped)){
			this.quality_string_encoded = this.rec.getAttribute("CQ").toString();
			this.quality_string_encoded = this.quality_string_encoded.
					substring(this.quality_string_encoded.length() - this.clipped_length, this.quality_string_encoded.length());
			this.quality_vector_decoded = QualityProcessing.decodeQualStringToVector(this.quality_string_encoded, phred_offset);
		}
	}
	
	public boolean isStartClipped(){
		return this.start_clipped;
	}
	
	public int getClippedLength(){
		return this.clipped_length;
	}
	
	public String getClippedQualityString(){
		return this.quality_string_encoded;
	}
	
	public Vector<Integer> getClippedQualityVector(){
		return this.quality_vector_decoded;
	}
	
	public int getClippedPosition(){
		if(this.start_clipped){
			return this.rec.getAlignmentStart();
		}else{
			return this.rec.getAlignmentEnd();
		}
	}
	
	public SAMRecord getSAMRecord(){
		return this.rec;
	}
	
	public String getPrimaryMobileMapping(){
		MobileSAMTag mobileTag = new MobileSAMTag();
		String mobileMapping = "";

		mobileTag.parse(this.rec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
		mobileMapping = mobileTag.getMobileCategoryNames().get(0);

		return mobileMapping;
	}
	
	public String getHomoPolymerMappingBasedOnMobileSAMTag(){
		MobileSAMTag mobileTag = new MobileSAMTag();


		if (this.rec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString() != null){
			mobileTag.parse(this.rec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			return mobileTag.getHomoPolymer();
		}

		return "";
	}
	
	public boolean hasCertainHomoPolymerMappingBasedOnMobileSAMTag(String polymer){
		MobileSAMTag mobileTag = new MobileSAMTag();

		if (this.rec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString() != null){
			mobileTag.parse(this.rec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			return (polymer.equalsIgnoreCase(mobileTag.getHomoPolymer()));
		}
		
		
		return false;	
	}
}
