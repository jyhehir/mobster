package org.umcn.me.samexternal;

import java.util.Arrays;
import java.util.Vector;


import net.sf.samtools.*;

/**
 * 
 * @author Djie
 *
 */
public class BAMReadFilter {

	//TODO: nr of mappings of read filter will not work on all SAMRecords
	//mapping tool has to be defined i.e. bwa or mosaik
	//and a subclass for that certain mapping tool needs to exist for NrMappingsSAMRecordHolder
	
	private int mapq = 0;
	private double average_qual = 0.0;
	private int nr_bases_with_min_qual = 0;
	private int min_qual = 0;
	private int max_edit_distance = 10000000;
	private int max_nr_mappings = 10000000;
	private int min_nr_mappings = 0;
	private String[] reference_names = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM", "1",
										"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"};
	private Vector<String> reference_regions = new Vector<String>(Arrays.asList(this.reference_names));
	private String mapping_tool = "";
	
	public BAMReadFilter(){
		
	}
	
	public BAMReadFilter(String mapping_tool){
		setMappingTool(mapping_tool);
	}
	
	
	
	public boolean passesFilter(SAMRecord record, boolean doNotTestForInReference){
		
		//Little bit dirty: using try catch for flow control

		try {
			NrMappingsSAMRecordHolder recHolder = NrMappingsSAMRecordHolderFactory.makeNrMappingsSAMRecordHolder(record, this.mapping_tool);
			return (recHolder.hasMAPQHigherThan(mapq) &&
					 recHolder.hasMeanQualHigherThan(average_qual) &&
					 recHolder.hasXbasesHigherThanQualY(nr_bases_with_min_qual, min_qual) &&
					 recHolder.hasLowerEditDistanceThan(max_edit_distance) &&
					 (doNotTestForInReference || recHolder.isInReferenceList(reference_regions)) &&
					 recHolder.getNumberOfMappings() <= this.max_nr_mappings &&
					 recHolder.getNumberOfMappings() >= this.min_nr_mappings);
		} catch (UnknownParamException e) {
			System.err.println("WARNING: can not filter on number of mappings!");
			SpecificSAMRecordHolder recHolder = new SpecificSAMRecordHolder(record);
			return (recHolder.hasMAPQHigherThan(mapq) &&
					 recHolder.hasMeanQualHigherThan(average_qual) &&
					 recHolder.hasXbasesHigherThanQualY(nr_bases_with_min_qual, min_qual) &&
					 recHolder.hasLowerEditDistanceThan(max_edit_distance) &&
					 (doNotTestForInReference || recHolder.isInReferenceList(reference_regions)));
		}

		
	}
	
	public Vector<String> regionsMatchingRecord(SAMRecord record, Vector<String> regions){
		Vector<String> matchingRegions = new Vector<String>();
		
		String ref = record.getReferenceName();
		SpecificSAMRecordHolder recHolder = new SpecificSAMRecordHolder(record);
		
		for (String region : regions){
			if (ref.equals(region)){
				matchingRegions.add(region);
			}else{
				if(region.contains(":")){
					String[] firstSplit = region.split(":");
					String refRegion = firstSplit[0];
					int startRegion = Integer.parseInt(firstSplit[1].split("-")[0]);
					int endRegion = Integer.parseInt(firstSplit[1].split("-")[1]);
					
					if(recHolder.isOnReferenceAndWithinBounds(refRegion, startRegion, endRegion)){
						matchingRegions.add(region);
					}
				}
			}
		}
		
		return matchingRegions;
	}
	
	public BAMReadFilter setMAPQFilter(int threshold){
		this.mapq = threshold;
		return this;
	}
	
	public BAMReadFilter setAverageQualFilter(double threshold){
		this.average_qual = threshold;
		return this;
	}
	
	public BAMReadFilter setXBasesToMinYQualFilter(int xBases, int yQual){
		this.nr_bases_with_min_qual = xBases;
		this.min_qual = yQual;
		return this;
	}
	
	public BAMReadFilter setMappingTool(String tool){
		this.mapping_tool = tool;
		return this;
	}
	
	public BAMReadFilter setNrOfMaxMappings(int max){
		this.max_nr_mappings = max;
		return this;
	}
	
	public BAMReadFilter setNrOfMinMappings(int min){
		this.min_nr_mappings = min;
		return this;
	}
	
	public BAMReadFilter setMaxEditDistance(int threshold){
		this.max_edit_distance = threshold;
		return this;
	}
	
	public BAMReadFilter setReferenceRegions(Vector<String> regions){
		this.reference_regions = regions;
		return this;
	}
	
	public String getNrOfMaxMappings(){
		return Integer.toString(this.max_nr_mappings);
	}
	
	public String getNrOfMinMappings(){
		return Integer.toString(this.min_nr_mappings);
	}
	
	public String getMaxEditDistance(){
		return Integer.toString(this.max_edit_distance);
	}
	
	public String getMinimumMAPQ(){
		return Integer.toString(this.mapq);
	}
	
	public String getMappingTool(){
		return this.mapping_tool;
	}
	
	
		
	
	
}
