package org.umcn.me.sam;


import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.umcn.me.samexternal.QualityProcessing;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.SampleBam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;



public class MateCluster<T extends SAMRecord> extends Vector<T> {
	
	private static final long serialVersionUID = 2309445895755132334L;
	
	private boolean assume_sorted;
	private SAMRecord record;
	private boolean split_read;
	private int multiple_hits;
	private int unique_hits;
	private int unmapped_hits;
	private Map<String, Integer> sample_count;
	private SampleBam sample_calling;
	
	public MateCluster(SAMFileHeader header, boolean split, boolean sorted, SampleBam sample){
		
		super();
		this.sample_calling = sample;
		this.assume_sorted = sorted;
		this.record = new SAMRecord(header);
		this.split_read = split;
		this.sample_count = new HashMap<String, Integer>();
		
	}
	
	public boolean add(T record){
		if (this.size() == 0){
			super.add(record);
			return true;
		}else{
			try {
				if (isValidAddition(this.firstElement(), record)){
					super.add(record);
					return true;
				}else{
					return false;
				}
			} catch (InvalidCategoryException e) {
				System.err.println("Could not add: " + record.getReadName() +
						"to cluster, because of invalid ME tag.");
				return false;
			}
		}
	}
	
	private void countCategories(){
		this.unique_hits = 0;
		this.multiple_hits = 0;
		this.unmapped_hits = 0;
		
		for(T record : this){
			String readName = record.getReadName();
			if (readName.startsWith(SAMDefinitions.UNIQUE_MULTIPLE_MAPPING)){
				this.multiple_hits += 1;
			}else if (readName.startsWith(SAMDefinitions.UNIQUE_UNIQUE_MAPPING)){
				this.unique_hits += 1;
			}else if (readName.startsWith(SAMDefinitions.UNIQUE_UNMAPPED_MAPPING)){
				this.unmapped_hits += 1;
			}
		}
	}
	
	private boolean isValidAddition(T originalRec, T newRec) throws InvalidCategoryException{
		if(!this.split_read){
			return (isSameStrand(originalRec, newRec) && isSameReference(originalRec, newRec) &&
			isSameMobileMapping(originalRec, newRec));
		}else{
			return (isSameReference(originalRec, newRec) &&	isSameMobileMapping(originalRec, newRec));
		}
	}
	
	private boolean isSameStrand(T originalRec, T newRec){
		return (originalRec.getReadNegativeStrandFlag() == newRec.getReadNegativeStrandFlag());
	}
	
	private boolean isSameReference(T originalRec, T newRec){
		return (originalRec.getReferenceName() == newRec.getReferenceName());
	}
	
	//TODO now this function returns only true when the first added mobile category
	//for two SAMRecords is the same (usually the best hit). Could make this more lenient
	//to i.e. also return true when the second mobile category of rec 1 is the same as the
	//first mobile category of rec 2.
	private boolean isSameMobileMapping(T originalRec, T newRec) throws InvalidCategoryException{
		MobileSAMTag originalTag = new MobileSAMTag();
		MobileSAMTag newTag = new MobileSAMTag();
		
		originalTag.parse(originalRec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
		newTag.parse(newRec.getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
		
		return (originalTag.getMobileCategoryNames().get(0).equals(newTag.getMobileCategoryNames().get(0)));
		
	}
	
	public boolean isWithinSearchArea(T record, int searchArea){
		
		boolean sameReference;
		boolean sameRange;
		
		if(this.size() == 0){
			return true;
		}else{
			sameReference = (record.getReferenceName() == this.lastElement().getReferenceName());
			sameRange = (record.getAlignmentStart() <= this.lastElement().getAlignmentStart() + searchArea);
			
			return (sameReference && sameRange);
			
			//Old return method in 0.1.6
			//return (record.getAlignmentStart() <= this.lastElement().getAlignmentStart() + searchArea);
		}
	}
	
	public int getClusterStart(){
		int start = 0;
		
		if (this.assume_sorted){
			start = this.firstElement().getAlignmentStart();
		}
		
		return start;
	}
	
	public int getClusterEnd(){
		int end = 0;
		int currentEnd = 0;
		
		for (SAMRecord rec : this){
			currentEnd = rec.getAlignmentEnd();
			if (currentEnd > end){
				end = currentEnd;
			}
		}
		
		return end;
	}
	
	private void countSamples(){
		String name;
		int old_count;
		
		for (SAMRecord rec : this){
			if (this.sample_calling.equals(SampleBam.MULTISAMPLE)){
				if(rec.getReadGroup() == null){
					System.err.println("[MateCluster] WARNING Multisample calling enabled but read does not contain @RG. Read will be counted as 'NoRG': " + rec.getReadName());
					name = "NoRG";
				}else{
					name = rec.getReadGroup().getSample();
				}
			}else if(this.sample_calling.equals(SampleBam.SINGLESAMPLE)){
				name = rec.getAttribute(MobileDefinitions.SAM_TAG_SAMPLENAME).toString();
			}else{
				name = rec.getAttribute(MobileDefinitions.SAM_TAG_SAMPLENAME).toString();
			}
			if (this.sample_count.containsKey(name)){
				old_count = this.sample_count.get(name);
				this.sample_count.put(name, old_count + 1);
			}else{
				this.sample_count.put(name, 1);
			}
		}
	}
	
	public int getClusterSize(){
		return getClusterEnd() - getClusterStart() + 1;
	}
	
	public void writeClusterToSAMWriter(SAMFileWriter writer, String name){
		StringBuilder cigar = new StringBuilder();
		StringBuilder referenceBuilder = new StringBuilder();
		MobileSAMTag meTag = new MobileSAMTag();
		String reference = this.firstElement().getReferenceName();
		

			meTag.parse(this.firstElement().getAttribute(MobileDefinitions.SAM_TAG_MOBILE).toString());
			
			if (!reference.startsWith("chr")){
				referenceBuilder.append("chr");
			}
			referenceBuilder.append(reference);
			
			int clusterSize = this.getClusterSize();
			cigar.append(clusterSize);
			cigar.append("M");

			record.setReadName(name);
			record.setFlags(0);
			record.setReferenceName(referenceBuilder.toString());
			record.setAlignmentStart(this.getClusterStart());
			record.setMappingQuality(255);
			record.setCigarString(cigar.toString());
			record.setMateReferenceName("*");
			record.setInferredInsertSize(0);
			record.setReadString(QualityProcessing.createNSequence(clusterSize));
			record.setBaseQualityString("*");
			
			if(!this.firstElement().getReadNegativeStrandFlag()){
				//positive strand mapping
				record.setFlags(0);
			}else{
				//negative strand mapping
				record.setFlags(16);
			}
			
			record.setAttribute(MobileDefinitions.SAM_TAG_CLUSTER_HITS, Integer.toString(this.size()));
			record.setAttribute(MobileDefinitions.SAM_TAG_CLUSTER_LENGTH, Integer.toString(clusterSize));
			record.setAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT, meTag.getMobileCategoryNames().get(0));
			record.setAttribute(MobileDefinitions.SAM_TAG_SPLIT_CLUSTER, Boolean.toString(this.split_read));
			
			this.countCategories();
			record.setAttribute(MobileDefinitions.SAM_TAG_UNIQUE_HITS, this.unique_hits);
			record.setAttribute(MobileDefinitions.SAM_TAG_MULTIPLE_HITS, this.multiple_hits);
			record.setAttribute(MobileDefinitions.SAM_TAG_UNMAPPED_HITS, this.unmapped_hits);
			
			this.countSamples();
			String sampleCount = this.sample_count.toString();
			record.setAttribute(MobileDefinitions.SAM_TAG_SAMPLECOUNT, sampleCount.substring(1, sampleCount.length() - 1));
			
			
			cigar.setLength(0);
			
			writer.addAlignment(record);

		
	}
	

}
