package org.umcn.me.util;

import net.sf.samtools.SAMRecord;

public class ReadName {
	
	public final String readName;
	public final String reference;
	public final int start;
	public final int end;
	
	public final String mateReference;
	public final int mateStart;
	public final int mateEnd;
	
	private ReadNameOption option;
	
	public ReadName(SAMRecord record, ReadNameOption option){
		
		this.option = option;
		
		this.readName = getReadName(record, this.option.prefixLength);
		this.reference = option.prefixReference + record.getReferenceName();
		this.start = record.getAlignmentStart();
		this.end = record.getAlignmentEnd();
		
		this.mateReference = option.prefixReference + record.getMateReferenceName();
		this.mateStart = record.getMateAlignmentStart();
		this.mateEnd = record.getMateAlignmentStart() + record.getReadLength() - 1;

	}
	
	public String getReadName(SAMRecord record, int removePrefix){
		String name = record.getReadName();
		return name.substring(removePrefix, name.length());
	}
	
	
	public String toString(){
		if (!this.option.addRegion){
			return this.readName;
		}else{
			return this.readName + "\t" + this.reference + "\t" + 
					this.start + "\t" + this.end + "\t" + this.mateReference + "\t" + this.mateStart + "\t" + this.mateEnd;
		}
	}
	
	
}
