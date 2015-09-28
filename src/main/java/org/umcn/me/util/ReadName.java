package org.umcn.me.util;

import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

import org.umcn.me.tabix.BlacklistAnnotation;
import org.umcn.me.tabix.RefGeneAnnotation;
import org.umcn.me.tabix.RepMaskAnnotation;
import org.umcn.me.tabix.SelfChainAnnotation;

public class ReadName {
	
	public final static String STANDARD_REFERENCE_PREFIX = "chr";
	
	public final String readName;
	public final String reference;
	public final int start;
	public final int end;
	public final boolean isMapped;
	
	/**
	 * absolute value of insertSize
	 */
	public final int insertSize;
	
	public final String mateReference;
	public final int mateStart;
	public final int mateEnd;
	public final boolean mateIsMapped;
	
	private List<RepMaskAnnotation> repmask = new ArrayList<RepMaskAnnotation>();
	private List<RefGeneAnnotation> refgene = new ArrayList<RefGeneAnnotation>();
	private List<BlacklistAnnotation> blacklist = new ArrayList<BlacklistAnnotation>();
	private List<SelfChainAnnotation> selfChain = new ArrayList<SelfChainAnnotation>();
	
	private List<RepMaskAnnotation> mateRepMask = new ArrayList<RepMaskAnnotation>();
	private List<RefGeneAnnotation> mateRefGene = new ArrayList<RefGeneAnnotation>();
	private List<BlacklistAnnotation> mateBlacklist = new ArrayList<BlacklistAnnotation>();
	
	
	
	private ReadNameOption option;
	
	public ReadName(SAMRecord record, ReadNameOption option){
		
		this.option = option;
		
		this.insertSize = Math.abs(record.getInferredInsertSize());
		
		this.readName = getReadName(record, this.option.prefixLength);
		this.reference = this.getReferenceName(record.getReferenceName(), option);
		this.start = record.getAlignmentStart();
		this.end = record.getAlignmentEnd();
		this.isMapped = ! record.getReadUnmappedFlag();
		
		this.mateReference = this.getReferenceName(record.getMateReferenceName(), option);
		this.mateStart = record.getMateAlignmentStart();
		this.mateEnd = record.getMateAlignmentStart() + record.getReadLength() - 1;
		this.mateIsMapped = ! record.getMateUnmappedFlag();

	}
	
	public String getReadName(SAMRecord record, int removePrefix){
		String name = record.getReadName();
		return name.substring(removePrefix, name.length());
	}
	
	public String getReferenceName(String reference, ReadNameOption option){
		String ref = option.prefixReference + reference;
		if (option.autoPrefixReference){
			if (! ref.startsWith(STANDARD_REFERENCE_PREFIX)){
				ref = STANDARD_REFERENCE_PREFIX + ref;
			}
		}
		return ref;
	}
	
	public String toString(){
		if (!this.option.addRegion){
			return this.readName;
		}else{
			return this.readName + "\t" + this.reference + "\t" + this.start + "\t" + this.end + "\t" + this.isMapped + "\t" +
					this.mateReference + "\t" + this.mateStart + "\t" + this.mateEnd + "\t" + this.mateIsMapped;
		}
	}
	
	public String toPositionString(){
		return this.reference + ":" + this.start + "-" + this.end;
	}
	
	public String mateToPositionString(){
		return this.mateReference + ":" + this.mateStart + "-" + this.mateEnd;
	}
	
	public List<RepMaskAnnotation> getRepMaskAnnotation(){
		return this.repmask;
	}
	
	public List<RefGeneAnnotation> getRefGeneAnnotation(){
		return this.refgene;
	}
	
	public void setRepMaskAnnotation(List<RepMaskAnnotation> repmaskAnnotations){
		this.repmask = repmaskAnnotations;
	}
	
	public void setRefGeneAnnotation(List<RefGeneAnnotation> refGeneAnnotations){
		this.refgene = refGeneAnnotations;
	}
	
	public List<RepMaskAnnotation> getMateRepMaskAnnotation(){
		return this.mateRepMask;
	}
	
	public List<RefGeneAnnotation> getMateRefGeneAnnotation(){
		return this.mateRefGene;
	}
	
	public void setMateRepMaskAnnotation(List<RepMaskAnnotation> repmaskAnnotations){
		this.mateRepMask = repmaskAnnotations;
	}
	
	public void setMateRefGeneAnnotation(List<RefGeneAnnotation> refGeneAnnotations){
		this.mateRefGene = refGeneAnnotations;
	}

	public List<BlacklistAnnotation> getBlacklist() {
		return blacklist;
	}

	public void setBlacklist(List<BlacklistAnnotation> blacklist) {
		this.blacklist = blacklist;
	}

	public List<BlacklistAnnotation> getMateBlacklist() {
		return mateBlacklist;
	}

	public void setMateBlacklist(List<BlacklistAnnotation> mateBlacklist) {
		this.mateBlacklist = mateBlacklist;
	}

	public List<SelfChainAnnotation> getMateSelfChain() {
		return this.selfChain;
	}

	public void setSelfChain(List<SelfChainAnnotation> selfChain) {
		this.selfChain = selfChain;
	}
	
	public SelfChainAnnotation returnSelfChainOverlappingMate(){
		if (this.selfChain == null){
			return null;
		}
		
		for (SelfChainAnnotation chain : this.selfChain){
			SimpleRegion chainRegion = chain.toRegion();
			if (this.mateIsMapped){
				SimpleRegion mateRegion = new SimpleRegion(this.mateReference, this.mateStart, this.mateEnd);
				if (chainRegion.hasOverlapInBP(mateRegion) > 0){
					return chain;
				}
			}
		}
		return null;
	}
	
	public boolean selfChainOverlapsMate(){
		if (this.returnSelfChainOverlappingMate() == null){
			return false;
		}
		return true;
	}
	
	public String overlappingSelfChain(){
		SelfChainAnnotation chain = this.returnSelfChainOverlappingMate();
		if (chain == null){
			return "";
		}
		return chain.toRegion().toString();
	}
	
	public double returnScoreOfOverlappingSelfChain(){
		SelfChainAnnotation chain = this.returnSelfChainOverlappingMate();
		if (chain == null){
			return 0;
		}
		return chain.normScore;
	}
	
	
}
