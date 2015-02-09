package org.umcn.me.samexternal;

import java.io.PrintWriter;

import org.umcn.me.parsers.FastaParserFast;

import net.sf.samtools.SAMRecord;

/**
 * @author Djie
 */

public abstract class NrMappingsSAMRecordHolder extends SpecificSAMRecordHolder implements MappingCountable {
		
		protected int min_clipping = 1;
		protected int max_clipping = 100;
		protected int min_avg_qual = 20;
	
		public NrMappingsSAMRecordHolder(SAMRecord samRecord){
			super(samRecord);
		}
		
		public NrMappingsSAMRecordHolder(SAMRecord samRecord, int minClipping, int maxClipping){
			this(samRecord);
			this.min_clipping = minClipping;
			this.max_clipping = maxClipping;
		}
		
		public void writeSAMRecordToFastQPrintWriter(PrintWriter outFq, boolean appendReadNumber,
				boolean includeSplit){
			StringBuilder readName = new StringBuilder(this.sam_record.getReadName());
			String sequence;
			String qualityString;			
			

			if(includeSplit && this.isSoftClippedMappedProperAndCertainSizeAndCertainQual()){
				sequence = this.getLargestClippedSeq();
				qualityString = this.getLargestClippedQual();
				
			}else{
				sequence = this.sam_record.getReadString();
				qualityString = this.sam_record.getBaseQualityString();
			}
			
			if(appendReadNumber && this.sam_record.getReadPairedFlag()){
				readName.append(getReadNumber());
			}
			
			FastaParserFast.writeFastQToStream(outFq, readName.toString(), sequence, qualityString);
		}
		
		public void setMinAvgQual(int qual){
			this.min_avg_qual = qual;
		}
		
		private String getReadNumber(){
			StringBuilder readSuffix = new StringBuilder(SAMDefinitions.READ_NUMBER_SEPERATOR);
			if (this.sam_record.getFirstOfPairFlag()){
				readSuffix.append("1");
				return readSuffix.toString();
			} else if (this.sam_record.getSecondOfPairFlag()){
				readSuffix.append("2");
				return readSuffix.toString();
			} else{
				return "";
			}
		}
		
		public boolean isSoftClippedAndMappedProper(){
			if (this.soft_clipping_left > 0 || this.soft_clipping_right > 0){
				if (this.isMappedUniquely() || this.sam_record.getProperPairFlag()){
					return true;
				}
			}
			
			return false;
		}
		
		public boolean isSoftClippedMappedProperAndCertainSizeAndCertainQual(){
			
			if (QualityProcessing.calculateAverage(QualityProcessing.decodeQualStringToVector(this.getLargestClippedQual(), 33))
					< min_avg_qual ){
				return false;
			}
			
			if (this.isMappedUniquely() || this.sam_record.getProperPairFlag()) {
				if (this.soft_clipping_left <= this.max_clipping
						&& this.soft_clipping_right >= this.min_clipping) {
					return true;
				} else if (this.soft_clipping_right <= this.max_clipping
						&& this.soft_clipping_left >= this.min_clipping) {
					return true;
				}
			}
			return false;
		}
}
