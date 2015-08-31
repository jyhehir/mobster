package org.umcn.me.sam;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.NrMappingsSAMRecordHolder;
import org.umcn.me.samexternal.NrMappingsSAMRecordHolderFactory;
import org.umcn.me.samexternal.SAMRecordHolderPair;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.samexternal.UnknownParamException;

import com.google.code.jyield.Generator;
import com.google.code.jyield.Yieldable;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class PotentialMobilePairIterator implements Generator<SAMRecordHolderPair<NrMappingsSAMRecordHolder>>{
	
	private SAMFileReader sam_reader = null;
	private String mapping_tool;
	Map<String, SAMRecord> potential_mobile_read_SAMRec_map = null;
	
	private int min_clipping_split = 35;
	private int max_clipping_split_other_side = 7;
	private int min_avg_qual = 0;
	
	private boolean include_split;
	private int min_anchor_mapq;
	private boolean skip_um_pairs = false; //unique -multiple pairs
	
	public static Logger logger = Logger.getLogger("PotentialMobilePairIterator");
	
	
	public PotentialMobilePairIterator(File samFile, String tool, boolean includeSplitReads){
		this.sam_reader = new SAMSilentReader(samFile);
		this.mapping_tool = tool;
		this.potential_mobile_read_SAMRec_map = new HashMap<String, SAMRecord>();
		this.include_split = includeSplitReads;
	}
	
	public PotentialMobilePairIterator(File samFile, String tool, boolean includeSplitReads,
			int minClipping, int maxClipping, int minQual, int minMapq){
		this(samFile, tool, includeSplitReads);
		this.min_clipping_split = minClipping;
		this.max_clipping_split_other_side = maxClipping;
		this.min_avg_qual = minQual;
		this.min_anchor_mapq = minMapq;
	}
	
	public SAMFileReader getSAMReader(){
		return this.sam_reader;
	}
	
	public void setSkippingOfUMPairs(boolean skip){
		this.skip_um_pairs = skip;
	}
	
	public void generate(Yieldable<SAMRecordHolderPair<NrMappingsSAMRecordHolder>> yieldable){
		
		try {
			for (SAMRecord rec : this.sam_reader){
				yieldPotentialMobileReadPair(rec, yieldable);				
				
			}
		} catch (UnknownParamException e) {
			logger.error("Unknown Mapping tool : " + e.getMessage());
			closeSamReader();
		} catch (IllegalSAMPairException e){
			logger.error(e.getMessage());
			closeSamReader();
		}
		closeSamReader();
	}

	private void yieldPotentialMobileReadPair(SAMRecord rec, Yieldable<SAMRecordHolderPair<NrMappingsSAMRecordHolder>> yieldable)
			throws UnknownParamException, IllegalSAMPairException {
		
		String readName;
		
		NrMappingsSAMRecordHolder potentialMobileRead1;
		NrMappingsSAMRecordHolder potentialMobileRead2;
		SAMRecordHolderPair<NrMappingsSAMRecordHolder> potentialReadPair;
		
		//Skip reads which are unmapped
		if((rec.getReadUnmappedFlag() && rec.getMateUnmappedFlag())){
			return;
		}
		
		//Skip single-end reads for now; will hog all the memory
		if (!rec.getReadPairedFlag()){
			return;
		}
		
		//In addition skip reads which are flagged as duplicate, are not primary or are supplemental
		//Note this will lead to memory hogging if from a pair one end is marked as duplicate and the other end isn't.
		//Although I think this should not happen.
		if (rec.getDuplicateReadFlag() || rec.getNotPrimaryAlignmentFlag() || rec.getSupplementaryAlignmentFlag()){
			return;
		}
		
		readName = rec.getReadName();
		
		if(!potential_mobile_read_SAMRec_map.containsKey(readName)){
			potential_mobile_read_SAMRec_map.put(readName, rec);
		}else{
			potentialMobileRead1 = NrMappingsSAMRecordHolderFactory.makeNrMappingsSAMRecordHolder(
					potential_mobile_read_SAMRec_map.get(readName), this.mapping_tool, this.min_clipping_split,
					this.max_clipping_split_other_side, this.min_anchor_mapq);
			potentialMobileRead2 = NrMappingsSAMRecordHolderFactory.makeNrMappingsSAMRecordHolder(
					rec, this.mapping_tool, this.min_clipping_split, this.max_clipping_split_other_side, this.min_anchor_mapq);
			
			if (this.min_avg_qual != 0){
				potentialMobileRead1.setMinAvgQual(this.min_avg_qual);
				potentialMobileRead2.setMinAvgQual(this.min_avg_qual);
			}
			
			potentialReadPair = new SAMRecordHolderPair<NrMappingsSAMRecordHolder>(potentialMobileRead1, potentialMobileRead2);
			
			if (!this.skip_um_pairs){
				if (potentialReadPair.isPotentialMobilePair(this.include_split)){
					yieldable.yield(new SAMRecordHolderPair<NrMappingsSAMRecordHolder>(potentialMobileRead1, potentialMobileRead2));
				}
			}else{
				if (potentialReadPair.isPotentialGRIPPair(this.include_split)){
					yieldable.yield(new SAMRecordHolderPair<NrMappingsSAMRecordHolder>(potentialMobileRead1, potentialMobileRead2));
				}
			}
					
			potential_mobile_read_SAMRec_map.remove(readName);
		}
	}
	
	private void closeSamReader(){
		this.sam_reader.close();
	}
	

}
