package org.umcn.me.util;

import net.sf.samtools.SAMRecord;

import java.util.List;

import org.umcn.me.tabix.RefGeneAnnotation;
import org.umcn.me.tabix.RepMaskAnnotation;

public class ReadName {
	
	public final String readName;
	public final String reference;
	public final int start;
	public final int end;
	public final boolean isMapped;
	
	public final String mateReference;
	public final int mateStart;
	public final int mateEnd;
	public final boolean mateIsMapped;
	
	private List<RepMaskAnnotation> repmask;
	private List<RefGeneAnnotation> refgene;
	
	private List<RepMaskAnnotation> mateRepMask;
	private List<RefGeneAnnotation> mateRefGene;
	
	private ReadNameOption option;
	
	public ReadName(SAMRecord record, ReadNameOption option){
		
		this.option = option;
		
		this.readName = getReadName(record, this.option.prefixLength);
		this.reference = option.prefixReference + record.getReferenceName();
		this.start = record.getAlignmentStart();
		this.end = record.getAlignmentEnd();
		this.isMapped = ! record.getReadUnmappedFlag();
		
		this.mateReference = option.prefixReference + record.getMateReferenceName();
		this.mateStart = record.getMateAlignmentStart();
		this.mateEnd = record.getMateAlignmentStart() + record.getReadLength() - 1;
		this.mateIsMapped = ! record.getMateUnmappedFlag();

	}
	
	public String getReadName(SAMRecord record, int removePrefix){
		String name = record.getReadName();
		return name.substring(removePrefix, name.length());
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
	
	
}
