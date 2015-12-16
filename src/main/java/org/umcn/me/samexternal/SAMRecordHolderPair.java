package org.umcn.me.samexternal;

import java.io.PrintWriter;
import java.util.Vector;



import net.sf.samtools.SAMFileWriter;

/**
 * @author Djie
 */

public class SAMRecordHolderPair<T extends NrMappingsSAMRecordHolder>{
	
	private NrMappingsSAMRecordHolder read1;
	private NrMappingsSAMRecordHolder read2;
	
	public SAMRecordHolderPair(T read1, T read2) throws IllegalSAMPairException{
		if(read1.getSAMRecord().getReadName().equals(read2.getSAMRecord().getReadName())){
			this.read1 = read1;
			this.read2 = read2;
		}else{
			throw new IllegalSAMPairException(read1.getSAMRecord().getReadName() + " is not a mate of " +
					read2.getSAMRecord().getReadName());
		}

	}
	
	public NrMappingsSAMRecordHolder getRead1(){
		return this.read1;
	}
	
	public NrMappingsSAMRecordHolder getRead2(){
		return this.read2;
	}
	
	public Vector<NrMappingsSAMRecordHolder> getReads(){
		Vector<NrMappingsSAMRecordHolder> reads = new Vector<NrMappingsSAMRecordHolder>();
		reads.add(this.read1);
		reads.add(this.read2);
		return reads;
	}
	
	public boolean hasAnchor(){
		if (this.read1.isMappedUniquely() || this.read2.isMappedUniquely()){
			return true;
		}
		return false;
	}
	
	public boolean isPotentialMobilePair(boolean includeSplit){
		if (this.isMappedUUDiscordantly()){
			return true;
		}else if (this.hasUnmappedRead() && this.hasAnchor()){
			return true;
			//TODO: test leaving out this.hasMultipleMappedRead() --> when definining multiple and unique
			//using MAPQ this might not return all inproperly paired pairs.
		} else if (!this.read1.getSAMRecord().getProperPairFlag() && this.hasMultipleMappedRead() && this.hasAnchor()){
			return true;
		} else if (includeSplit && this.hasSplitReadOfCertainSize()){
			return true;
		}
		return false;
	}
	
	/**
	 * Unlike isPotentialMobilePair, does not check for UM pairs
	 * @param includeSplit
	 * @return
	 */
	public boolean isPotentialGRIPPair(boolean includeSplit){
		if (this.isMappedUUDiscordantly()){
			return true;
		}else if (this.hasUnmappedRead() && this.hasAnchor()){
			return true;
			//TODO: test leaving out this.hasMultipleMappedRead() --> when definining multiple and unique
			//using MAPQ this might not return all inproperly paired pairs.
		} else if (includeSplit && this.hasSplitReadOfCertainSize()){
			return true;
		}
		return false;
	}
	
	public boolean hasSplitReadOfCertainSize(){
		return (this.read1.isSoftClippedMappedProperAndCertainSizeAndCertainQual() ||
				this.read2.isSoftClippedMappedProperAndCertainSizeAndCertainQual());
				
	}
	
	public boolean hasMultipleMappedRead(){
		return (this.read1.isMappedMultiple() || this.read2.isMappedMultiple());
	}
	
	public boolean hasUnmappedRead(){
		return (this.read1.getSAMRecord().getReadUnmappedFlag() || this.read2.getSAMRecord().getReadUnmappedFlag());
	}
	
	public boolean isMappedUUDiscordantly(){
		return (this.read1.isMappedUniquely() && this.read2.isMappedUniquely() &&
				!this.read1.getSAMRecord().getProperPairFlag() && !this.read2.getSAMRecord().getProperPairFlag());
	}
	
	public boolean isMappedUniqueMultipleDiscordantly(){
		if (this.hasAnchor() && !this.read1.getSAMRecord().getProperPairFlag()){
			if (this.read1.isMappedMultiple() || this.read2.isMappedMultiple()){
				return true;
			}
		}
		return false;
	}
	
	public boolean isMappedUniqueUnmapped(){
		if (this.hasAnchor()){
			if (this.read1.getSAMRecord().getReadUnmappedFlag() || this.read2.getSAMRecord().getReadUnmappedFlag()){
				return true;
			}
		}
		return false;
	}
	
	public Vector<NrMappingsSAMRecordHolder> getPotentialMobileReads(boolean includeSplitReads){
		Vector<NrMappingsSAMRecordHolder> potentialMobileReads = new Vector<NrMappingsSAMRecordHolder>();
		
		
		//split reads have priority
		if(this.hasSplitReadOfCertainSize() && includeSplitReads){
			if(this.read1.isSoftClippedMappedProperAndCertainSizeAndCertainQual()){
				potentialMobileReads.add(read1);
			}else if (this.read2.isSoftClippedMappedProperAndCertainSizeAndCertainQual()){
				potentialMobileReads.add(read2);
			}
		}else if (isMappedUUDiscordantly()){
			potentialMobileReads.add(read1);
			potentialMobileReads.add(read2);
		}else if(this.isMappedUniqueMultipleDiscordantly()){
			if(this.read1.isMappedMultiple()){
				potentialMobileReads.add(read1);
			}else{
				potentialMobileReads.add(read2);
			}
		}else if (this.isMappedUniqueUnmapped()){
			if(this.read1.getSAMRecord().getReadUnmappedFlag()){
				potentialMobileReads.add(read1);
			}else{
				potentialMobileReads.add(read2);
			}
		}
		return potentialMobileReads;	
	}
	
	public boolean writePotentialMobileReadsToFastQ(PrintWriter out, boolean includeSplitReads){
		
		Vector<NrMappingsSAMRecordHolder> potentialMobileReads = this.getPotentialMobileReads(includeSplitReads);
		
		for (NrMappingsSAMRecordHolder read : potentialMobileReads){
			read.writeSAMRecordToFastQPrintWriter(out, true, includeSplitReads);
		}
		
		if (potentialMobileReads.size() > 0){
			return true;
		}else{
			return false;
		}
			
	}
	
	public void writeMobileReadsToFastQAndPotentialPairsToSAM(PrintWriter fq,
			SAMFileWriter samWriter, boolean includeSplitReads, boolean prefix, String rgPrefix, String readSuffix){
		
		if (prefix){
			prefixReadNames(includeSplitReads);
		}
		
		//Hack below for William Brandler to handle identical readnames within same BAM file from multiple samples
		if (! "".equals(readSuffix)){
			this.suffixReadNames(readSuffix);
		}
		
		if (this.writePotentialMobileReadsToFastQ(fq, includeSplitReads)){
			this.writeToSAMFileWriter(samWriter, rgPrefix);
		}
		
	}

	private void prefixReadNames(boolean includeSplitReads) {
		
		StringBuilder readName1 = new StringBuilder();
		StringBuilder readName2 = new StringBuilder();
		String originalRead1 = this.read1.getSAMRecord().getReadName();
		String originalRead2 = this.read2.getSAMRecord().getReadName();
		
		//Split reads have priority
		if (this.hasSplitReadOfCertainSize() && includeSplitReads){
			Vector<NrMappingsSAMRecordHolder> potentialSplit = this.getPotentialMobileReads(true);
			readName1.append(SAMDefinitions.SPLIT_MAPPING);
			readName2.append(SAMDefinitions.SPLIT_MAPPING);
			if(potentialSplit.firstElement().isLargetsClippedSeqLeft()){
				readName1.append(SAMDefinitions.LEFT_CLIPPED);
				readName2.append(SAMDefinitions.LEFT_CLIPPED);
			}else{
				readName1.append(SAMDefinitions.RIGHT_CLIPPED);
				readName2.append(SAMDefinitions.RIGHT_CLIPPED);
			}
		}else if(this.isMappedUniqueMultipleDiscordantly()){
			readName1.append(SAMDefinitions.UNIQUE_MULTIPLE_MAPPING);
			readName2.append(SAMDefinitions.UNIQUE_MULTIPLE_MAPPING);
		}else if(this.isMappedUniqueUnmapped()){
			readName1.append(SAMDefinitions.UNIQUE_UNMAPPED_MAPPING);
			readName2.append(SAMDefinitions.UNIQUE_UNMAPPED_MAPPING);
		}else if (this.isMappedUUDiscordantly()){
			readName1.append(SAMDefinitions.UNIQUE_UNIQUE_MAPPING);
			readName2.append(SAMDefinitions.UNIQUE_UNIQUE_MAPPING);
		}
		readName1.append(originalRead1);
		readName2.append(originalRead2);
		this.read1.getSAMRecord().setReadName(readName1.toString());
		this.read2.getSAMRecord().setReadName(readName2.toString());
	}
	
	private void suffixReadNames(String suffix){
		
		StringBuilder readName1 = new StringBuilder(this.read1.getSAMRecord().getReadName());
		StringBuilder readName2 = new StringBuilder(this.read2.getSAMRecord().getReadName());
		
		readName1.append(suffix);
		readName2.append(suffix);
		
		this.read1.getSAMRecord().setReadName(readName1.toString());
		this.read2.getSAMRecord().setReadName(readName2.toString());

	}
	
	public void writeToSAMFileWriter(SAMFileWriter samWriter, String rgPrefix){

		
		
		if (rgPrefix != null && ! rgPrefix.equals("")){
			//code for read1
			if (this.read1.getSAMRecord().getAttribute("RG") == null){
				this.read1.getSAMRecord().setAttribute("RG", rgPrefix);
			} else{
				StringBuilder readGroup1 = new StringBuilder(rgPrefix);
				readGroup1.append((String) this.read1.getSAMRecord().getAttribute("RG"));
				this.read1.getSAMRecord().setAttribute("RG", readGroup1.toString());
			}
			
			//code for read2
			if (this.read2.getSAMRecord().getAttribute("RG") == null){
				this.read2.getSAMRecord().setAttribute("RG", rgPrefix);
			}else{
				StringBuilder readGroup2 = new StringBuilder(rgPrefix);
				readGroup2.append((String) this.read2.getSAMRecord().getAttribute("RG"));
				this.read2.getSAMRecord().setAttribute("RG", readGroup2.toString());
			}
			
		}
		
		samWriter.addAlignment(this.read1.getSAMRecord());
		samWriter.addAlignment(this.read2.getSAMRecord());
	}


}
