package org.umcn.me.util;

import java.math.BigDecimal;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

import org.umcn.gen.sam.QualityProcessing;
import org.umcn.me.sam.InvalidCategoryException;
import org.umcn.me.sam.MobileSAMTag;
import org.umcn.me.splitread.ClippedRead;

/**
 * Vector Containing ClippedReads or subclasses of ClippedReads.
 * 
 * This class contains convenience methods to extract annotation information
 * about these clipped reads. For instance to get the average base
 * PHRED score for all clipped ends in the vector.
 * 
 * NOTE: All getter methods should speak for themselves and are therefore
 * not commented
 * 
 * 
 * @author Djie
 *
 * @param <T>
 */
public class ClippedReadSet<T extends ClippedRead> extends Vector<T>{

	private static final long serialVersionUID = -1303033868736365325L;
	
	public ClippedReadSet<ClippedRead> getClippedReadsAtOnePosition(boolean left){
		ClippedReadSet<ClippedRead> subset = new ClippedReadSet<ClippedRead>();
		
		for (ClippedRead c : this){
			if(left && c.isStartClipped()){
				subset.add(c);
			}else if(!left && !c.isStartClipped()){
				subset.add(c);
			}
		}
		return subset;
	}
	
	public Vector<Integer> getClippedLengths(){
		Vector<Integer> lengths = new Vector<Integer>();
		
		for (ClippedRead c : this){
			lengths.add(c.getClippedLength());
		}
		return lengths;
	}
	
	public Vector<Integer> getClippedLengthsLeft(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		return subset.getClippedLengths();
	}
	
	public Vector<Integer> getClippedLengthsRight(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		return subset.getClippedLengths();
	}
	
	public double getAverageClippedLength(){
		return QualityProcessing.calculateAverage(this.getClippedLengths());
	}
	
	public double getAverageClippedLengthLeft(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		if (subset.size() > 0){
			return QualityProcessing.calculateAverage(subset.getClippedLengths());
		}
		return -1;
	}
	
	public double getAverageClippedLengthRight(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		if (subset.size() > 0){
			return QualityProcessing.calculateAverage(subset.getClippedLengths());
		}
		return -1;
	}
	
	public double getAverageBaseQualAll(){
		double sum = 0.0;
		double avg;
		for(ClippedRead c : this){
			sum += QualityProcessing.calculateAverage(c.getClippedQualityVector());
		}
		
		avg = sum / this.size();
		
		return QualityProcessing.round(avg, 2, BigDecimal.ROUND_HALF_UP);
	}
	
	public double getAverageBaseQualLeftClipped(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		if (subset.size() > 0){
			return subset.getAverageBaseQualAll();
		}else{
			return -1;
		}
	}
	
	public double getAverageBaseQualRightClipped(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		if (subset.size() > 0){
			return subset.getAverageBaseQualAll();
		}else{
			return -1;
		}
	}
	
	public Vector<Integer> getClippedEnds(){
		Vector<Integer> clippedEnds = new Vector<Integer>();
		for(ClippedRead c : this){
			clippedEnds.add(c.getClippedPosition());
		}
		return clippedEnds;
	}
	
	public int getNumberLeftClippedReads(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		return subset.size();
	}
	
	public int getNumberRightClippedReads(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		return subset.size();
	}
	
	public int getMaxDistanceLeftClippedReads(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		if(subset.size() > 0){
			Vector<Integer> ends = subset.getClippedEnds();
			Collections.sort(ends);
			return ends.lastElement() - ends.firstElement();
		}
		return -1;
	}
	
	public int getMaxDistanceRightClippedReads(){
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		if(subset.size() > 0){
			Vector<Integer> ends = subset.getClippedEnds();
			Collections.sort(ends);
			return ends.lastElement() - ends.firstElement();
		}
		return -1;
	}
	
	public double getFractionOfLeftClippedReadsWithSamePos(){
		int max = 0;
		int temp;
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(true);
		if (subset.size() > 0){
			Vector<Integer> estimates = subset.getClippedEnds();
			for (int estimate : estimates){
				temp = Collections.frequency(estimates, estimate);
				if(temp > max){
					max = temp;
				}
			}
			return QualityProcessing.round((double) max / subset.size(),
					2, BigDecimal.ROUND_HALF_UP);
		}
		return 0.0;
	}
	
	public double getFractionOfRightClippedReadsWithSamePos(){
		int max = 0;
		int temp;
		ClippedReadSet<ClippedRead> subset = this.getClippedReadsAtOnePosition(false);
		if (subset.size() > 0){
			Vector<Integer> estimates = subset.getClippedEnds();
			for (int estimate : estimates){
				temp = Collections.frequency(estimates, estimate);
				if(temp > max){
					max = temp;
				}
			}
			return QualityProcessing.round((double) max / subset.size(),
					2, BigDecimal.ROUND_HALF_UP);
		}
		return 0.0;
	}
	
	public int getNrOfAluPredictions() throws InvalidCategoryException{
		int sum = 0;
		for(ClippedRead c : this){
			MobileSAMTag mobileTag = new MobileSAMTag();
			mobileTag.parse(c.getSAMRecord().getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			if ("ALU".equalsIgnoreCase(mobileTag.getMobileCategoryNames().get(0))){
				sum++;
			}
		}
		return sum;

	}
	
	public int getNrOfSVAPredictions() throws InvalidCategoryException{
		int sum = 0;
		for(ClippedRead c : this){
			MobileSAMTag mobileTag = new MobileSAMTag();
			mobileTag.parse(c.getSAMRecord().getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			if ("SVA".equalsIgnoreCase(mobileTag.getMobileCategoryNames().get(0))){
				sum++;
			}
		}
		return sum;
	}
	
	public int getNrOfL1Predictions() throws InvalidCategoryException{
		int sum = 0;
		for(ClippedRead c : this){
			MobileSAMTag mobileTag = new MobileSAMTag();
			mobileTag.parse(c.getSAMRecord().getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			if ("L1".equalsIgnoreCase(mobileTag.getMobileCategoryNames().get(0))){
				sum++;
			}
		}
		return sum;
	}
	
	public int getMedianOfLeftClippedEnds(){
		
		Vector<Integer> clippedEnds = this.getLeftClippedEnds();
		
		if(clippedEnds.size() != 0){
			Collections.sort(clippedEnds);
			return QualityProcessing.calculateMedian(clippedEnds);
		}
		
		return 0;		
	}
	
	public int getMedianOfRightClippedEnds(){
		
		Vector<Integer> clippedEnds = this.getRightClippedEnds();
		
		if(clippedEnds.size() != 0){
			Collections.sort(clippedEnds);
			return QualityProcessing.calculateMedian(clippedEnds);
		}
		
		return 0;		
	}
	
	public int getLowestClippedPosition(){
		Vector<Integer> clippedEnds = this.getClippedEnds();
		if (clippedEnds.size() > 0){
			Collections.sort(clippedEnds);
			return clippedEnds.firstElement();
		}
		return 0;
	}
	
	public int getHighestClippedPosition(){
		Vector<Integer> clippedEnds = this.getClippedEnds();
		if (clippedEnds.size() > 0){
			Collections.sort(clippedEnds);
			return clippedEnds.lastElement();
		}
		return 0;
	}
	
	public Vector<Integer> getLeftClippedEnds(){
		ClippedReadSet<ClippedRead> clippedReadsLeft = this.getClippedReadsAtOnePosition(true);
		Vector<Integer> clippedEnds = clippedReadsLeft.getClippedEnds();
		return clippedEnds;
	}
	
	public Vector<Integer> getRightClippedEnds(){
		ClippedReadSet<ClippedRead> clippedReadsLeft = this.getClippedReadsAtOnePosition(false);
		Vector<Integer> clippedEnds = clippedReadsLeft.getClippedEnds();
		return clippedEnds;
	}
	
	public int getDistanceBetweenLowestAndHighestClippedPosition(){
		return this.getHighestClippedPosition() - this.getLowestClippedPosition() + 1;
	}
	
	//lenient meaning: same primary mobile mapping is only required for
	//one end of the insertion not on both ends
	public boolean isSamePrimaryMobileMapping(ClippedRead read, boolean lenient){
		
		String primaryMobileMappingOfReadToBeAdded;
		String primaryMobileMappingOfFirstReadInSet;
		ClippedReadSet<ClippedRead> clippedReads;
		
		primaryMobileMappingOfReadToBeAdded = read.getPrimaryMobileMapping();
		
		//If query does not have primary mobile mapping, because of poly A/T mapping, return false
		if (primaryMobileMappingOfReadToBeAdded.equals("")){
			return false;
		}
		
		//If query has a primary mobile mapping which also occurs in the existing clippedreadset return true
		if (this.getMobileMapping().equals("") || this.getMobileMapping().contains(primaryMobileMappingOfReadToBeAdded)){
			return true;
		}
		
		
		if(this.size() > 0){
			
			if(!lenient){
				primaryMobileMappingOfFirstReadInSet = this.firstElement().getPrimaryMobileMapping();
			}else {
				clippedReads = this.getClippedReadsAtOnePosition(read.isStartClipped());
				if (clippedReads.size() > 0){
					primaryMobileMappingOfFirstReadInSet = clippedReads.firstElement().getPrimaryMobileMapping();
				}else{
					return true;
				}
			}			
			
			if(!primaryMobileMappingOfReadToBeAdded.equals(primaryMobileMappingOfFirstReadInSet)){
				return false;
			}
		}
		
		//If no read has been added to ClippedReadSet, return true;
		return true;
	}
	
	public boolean hasSameHomoPolymerMapping(ClippedRead read){
		
		if(read.getHomoPolymerMappingBasedOnMobileSAMTag().equals("")){
			return false;
		}
		
		if (this.size() > 0){
			ClippedReadSet<ClippedRead> clippedReads = this.getClippedReadsAtOnePosition(read.isStartClipped());
			if (clippedReads.size() > 0){
				String originalPolymerMapping = clippedReads.firstElement().getHomoPolymerMappingBasedOnMobileSAMTag();
				if (!read.hasCertainHomoPolymerMappingBasedOnMobileSAMTag(originalPolymerMapping)){
					return false;
				}
			}
		}
		
		return true;
	}
	
	public int getNrOfLeftClippedPolyAMapping(){
		ClippedReadSet<ClippedRead> clippedReads = this.getClippedReadsAtOnePosition(true);
		int nrOfMappings = 0;
		
		for (ClippedRead read : clippedReads){
			if (read.hasCertainHomoPolymerMappingBasedOnMobileSAMTag("polyA")){
				nrOfMappings++;
			}
		}
		return nrOfMappings;
	}
	
	public int getNrOfRightClippedPolyAMapping(){
		ClippedReadSet<ClippedRead> clippedReads = this.getClippedReadsAtOnePosition(false);
		int nrOfMappings = 0;
		
		for (ClippedRead read : clippedReads){
			if (read.hasCertainHomoPolymerMappingBasedOnMobileSAMTag("polyA")){
				nrOfMappings++;
			}
		}
		return nrOfMappings;
	}
	
	public int getNrOfLeftClippedPolyTMapping(){
		ClippedReadSet<ClippedRead> clippedReads = this.getClippedReadsAtOnePosition(true);
		int nrOfMappings = 0;
		
		for (ClippedRead read : clippedReads){
			if (read.hasCertainHomoPolymerMappingBasedOnMobileSAMTag("polyT")){
				nrOfMappings++;
			}
		}
		return nrOfMappings;
	}
	
	public int getNrOfRightClippedPolyTMapping(){
		ClippedReadSet<ClippedRead> clippedReads = this.getClippedReadsAtOnePosition(false);
		int nrOfMappings = 0;
		
		for (ClippedRead read : clippedReads){
			if (read.hasCertainHomoPolymerMappingBasedOnMobileSAMTag("polyT")){
				nrOfMappings++;
			}
		}
		return nrOfMappings;
	}
	
	public String getMobileMapping(){
		StringBuilder mobileFamily = new StringBuilder();
		String currentMobile;
		Set<String> mobileSet = new HashSet<String>();
		
		for (ClippedRead read : this){
			currentMobile = read.getPrimaryMobileMapping();
			if ("".equals(currentMobile)){
				continue;
			}
			if (!mobileSet.contains(currentMobile)){
				mobileSet.add(read.getPrimaryMobileMapping());
				mobileFamily.append(currentMobile);
				mobileFamily.append(";");
			}
		}
		
		if(mobileFamily.length() == 0){
			return "";
		}
		return mobileFamily.substring(0, mobileFamily.length() - 1);
		
	}
	
	
	public boolean addIfSameSecondaryMobileMapping(ClippedRead read){
		//TODO
		return false;
	}

	public boolean isRightClippedReadWithinRegion(ClippedRead read, int region){
		int clippedPos;
		int clippedPos2;
		ClippedReadSet<ClippedRead> rightClippedReads = this.getClippedReadsAtOnePosition(false);
		
		if(rightClippedReads.size() == 0){
			return true;
		}
		
		if (!read.isStartClipped()){
			clippedPos = read.getClippedPosition();
			clippedPos2 = rightClippedReads.firstElement().getClippedPosition();
			if (clippedPos >= clippedPos2 - region && clippedPos <= clippedPos2 + region){
				return true;
			}
		}
		return false;
	}
	
	public boolean isLeftClippedReadWithinRegion(ClippedRead read, int region){
		int clippedPos;
		int clippedPos2;
		ClippedReadSet<ClippedRead> leftClippedReads = this.getClippedReadsAtOnePosition(true);
		
		if(leftClippedReads.size() == 0){
			return true;
		}
		
		if(read.isStartClipped()){
			clippedPos = read.getClippedPosition();
			clippedPos2 = leftClippedReads.firstElement().getClippedPosition();
			if (clippedPos >= clippedPos2 - region && clippedPos <= clippedPos2 + region){
				return true;
			}
		}
		return false;
	}
	
	public boolean isOnSameReference(ClippedRead read){
		String oldReference;
		String newReference = read.getSAMRecord().getReferenceName();
		
		if(this.size() == 0){
			return true;
		}else{
			oldReference = this.firstElement().getSAMRecord().getReferenceName();
			if(newReference.equals(oldReference)){
				return true;
			}
		}
		return false;
	}
	
	public String getSampleCounts(){
		
		String name;
		Map<String, Integer> sampleCountMap = new HashMap<String, Integer>();
		int old_count;
		
		for (ClippedRead read : this) {
			name = read.getSAMRecord().getAttribute(MobileDefinitions.SAM_TAG_SAMPLENAME).toString();
			if (sampleCountMap.containsKey(name)) {
				old_count = sampleCountMap.get(name);
				sampleCountMap.put(name, old_count + 1);
			} else {
				sampleCountMap.put(name, 1);
			}
		}
		
		String sampleCountMapString = sampleCountMap.toString();
		
		return sampleCountMapString.substring(1, sampleCountMapString.length() - 1);
		
	}
	
	public void writeClippedReadSetToBAM(SAMFileWriter writer){
		
		StringBuilder cigar = new StringBuilder();
		SAMRecord record = new SAMRecord(writer.getFileHeader());
		int clusterSize = this.getDistanceBetweenLowestAndHighestClippedPosition();
		
		cigar.append(clusterSize);
		cigar.append("M");
		
		if (this.size()!= 0){
			record.setReadName(this.firstElement().getSAMRecord().getReadName());
			record.setFlags(0); //+ or - mapping does not matter in this case
			record.setReferenceName(this.firstElement().getSAMRecord().getReferenceName());
			record.setAlignmentStart(this.getLowestClippedPosition());
			record.setMappingQuality(255);
			record.setCigarString(cigar.toString());
			record.setMateReferenceName("*");
			record.setInferredInsertSize(0);
			record.setReadString(QualityProcessing.createNSequence(clusterSize));
			record.setBaseQualityString("*");
			setCustomSAMAttributes(record);
		}
		writer.addAlignment(record);
	}
	
	private void setCustomSAMAttributes(SAMRecord record){
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_CLIPPED_AVG_QUAL, Double.toString(this.getAverageBaseQualAll()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_AVG_CLIPPED_LEN, Double.toString(this.getAverageClippedLength()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_ENDS, CollectionUtil.toString(this.getLeftClippedEnds()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_FRAC_SAME_DISTANCE, Double.toString(this.getFractionOfLeftClippedReadsWithSamePos()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_HITS, this.getClippedReadsAtOnePosition(true).size());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_MAX_DISTANCE, this.getMaxDistanceLeftClippedReads());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_MEDIAN_END, this.getMedianOfLeftClippedEnds());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_POLYA, this.getNrOfLeftClippedPolyAMapping());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_POLYT, this.getNrOfLeftClippedPolyTMapping());
		
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_ENDS, CollectionUtil.toString(this.getRightClippedEnds()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_FRAC_SAME_DISTANCE, Double.toString(this.getFractionOfRightClippedReadsWithSamePos()));
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_HITS, this.getClippedReadsAtOnePosition(false).size());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_MAX_DISTANCE, this.getMaxDistanceRightClippedReads());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_MEDIAN_END, this.getMedianOfRightClippedEnds());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_POLYA, this.getNrOfRightClippedPolyAMapping());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_POLYT, this.getNrOfRightClippedPolyTMapping());
		
		record.setAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT, this.getMobileMapping());
		record.setAttribute(MobileDefinitions.SAM_TAG_SAMPLECOUNT, this.getSampleCounts());
		record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_CLUSTER, "true");
	}
	
	
	
	
}
