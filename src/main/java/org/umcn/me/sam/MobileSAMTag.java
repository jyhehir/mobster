package org.umcn.me.sam;

import java.util.Arrays;
import java.util.Vector;

import org.umcn.gen.sam.NrMappingsSAMRecordHolder;
import org.umcn.gen.sam.NrMappingsSAMRecordHolderFactory;
import org.umcn.gen.sam.UnknownParamException;
import org.umcn.me.util.MobileDefinitions;

import net.sf.samtools.SAMRecord;

public class MobileSAMTag {

	private String mobile_element = "";
	private String cigar = "";
	private String number_hits = "";
	private String homopolymer = "";
	private Vector<String> mobile_categories = null;
	private Vector<String> mobile_sam_tag = null;
	
	private final static String DELIMITER = ";";
	private final static String CATEGORY_DELIMITER = ",";
	
	private final int MOBILE_ELEMENT = 0;
	private final int CIGAR = 1;
	private final int MOBILE_CATEGORIES = 2;
	private final int MOBILE_HITS = 3;
	private final int HOMOPOLYMER = 4;
	
	
	public MobileSAMTag(){
		this.mobile_categories = new Vector<String>();
		this.mobile_sam_tag = new Vector<String>();
	}
		
	
	public void build(SAMRecord rec, String mapper) throws InvalidCategoryException, UnknownParamException{
		String reference = rec.getReferenceName();
		boolean unmapped = rec.getReadUnmappedFlag();
		
		for (int i = 0; i <= HOMOPOLYMER; i++){
			switch(i){
				case MOBILE_ELEMENT:
					if (!unmapped && "".equals(this.homopolymer)) setMobileElementName(reference);
					this.mobile_sam_tag.add(mobile_element);
					break;
				case CIGAR:
					if (!unmapped && "".equals(this.homopolymer)) setMobileCigar(rec.getCigarString());
					this.mobile_sam_tag.add(cigar);
					break;
				case MOBILE_CATEGORIES:
					if (!unmapped && "".equals(this.homopolymer)){
						addMobileCategoryByReference(reference);
					}else{
						this.mobile_categories.clear();
						this.mobile_categories.add("");
					}
					this.mobile_sam_tag.add(mobileCategoriesToString());
					break;
				case MOBILE_HITS:
					setNrOfMappings(rec, mapper);
					this.mobile_sam_tag.add(number_hits);
					break;
				case HOMOPOLYMER:
					this.mobile_sam_tag.add(this.homopolymer);
					break;
			}
		}
	}
	
	public void parse(String MEValueAttribute) throws InvalidCategoryException{
		String[] meValues = MEValueAttribute.split(DELIMITER, -1);
		
//		System.err.println(MEValueAttribute);
//		System.err.println(meValues.toString());
//		System.err.println(meValues[MOBILE_ELEMENT]);
//		System.err.println(meValues[CIGAR]);
//		System.err.println(meValues[MOBILE_CATEGORIES]);
//		System.err.println(meValues[MOBILE_HITS]);
//		System.err.println(meValues[HOMOPOLYMER]);
		
		this.setMobileElementName(meValues[MOBILE_ELEMENT]);
		this.setMobileCigar(meValues[CIGAR]);
		this.parseMobileCategories(meValues[MOBILE_CATEGORIES].split(CATEGORY_DELIMITER));
		this.setNrOfMappings(Integer.parseInt(meValues[MOBILE_HITS]));
		this.setHomoPolymer(meValues[HOMOPOLYMER]);
		
	}
	
	public void parseMobileCategories(String[] categories) throws InvalidCategoryException{
		for (String category : categories){
			this.addMobileCategory(category);
		}
	}
	
	public String mobileCategoriesToString(){
		StringBuilder sb = new StringBuilder();
		
		for(String category : this.mobile_categories){
			sb.append(category);
			sb.append(CATEGORY_DELIMITER);
		}
		
		return sb.substring(0, sb.length() - 1); //remove trailing category delimiter
	}


	public String toString(){
		StringBuilder mobileSAMValue = new StringBuilder();
		
		for (int i = 0; i < this.mobile_sam_tag.size(); i++){
			if (i == MOBILE_CATEGORIES){
				mobileSAMValue.append(this.mobileCategoriesToString());
			}else{
				mobileSAMValue.append(this.mobile_sam_tag.get(i));	
			}
			mobileSAMValue.append(DELIMITER);
		}
		
		return mobileSAMValue.toString();
	}


	public boolean addMobileCategory(String category) throws InvalidCategoryException{
		
		//The exception case where category == "", when read maps to polyA for instance
		if ("".equals(category)){
			if (this.mobile_categories.contains("")){
				return false;
			}else{
				return this.mobile_categories.add(category);
			}
		}
		
		for (String validCategory : MobileDefinitions.MOBILE_CATEGORIES){
			if (this.mobile_categories.contains(validCategory)){
				return false;
			}else if (validCategory.equalsIgnoreCase(category)){
				return this.mobile_categories.add(validCategory);
			}
		}
		throw new InvalidCategoryException("Wrong category: " + category + ". Must be one of: " + 
				Arrays.toString(MobileDefinitions.MOBILE_CATEGORIES));
	}


	public boolean addMobileCategoryByReference(String reference) throws InvalidCategoryException{
		
		for (String validCategory : MobileDefinitions.MOBILE_CATEGORIES){
			if (this.mobile_categories.contains(validCategory)){
				return false;
			}else if(reference.toUpperCase().startsWith(validCategory.toUpperCase())){
				return this.mobile_categories.add(validCategory);
			}
		}
		throw new InvalidCategoryException(reference + "can not be categorized into one of these" +
				" mobile categories: " + Arrays.toString(MobileDefinitions.MOBILE_CATEGORIES));
	}
	
	public String getHomoPolymer(){
		return this.homopolymer;
	}
	
	public String getMobileCigar(){
		return this.cigar;
	}
	
	public String getMobileElementName(){
		return this.mobile_element;
	}
	
	public Vector<String> getMobileCategoryNames(){
		return this.mobile_categories;
	}
	
	public int getNrOfMappings(){
		return Integer.parseInt(this.number_hits);
	}

	private void setMobileCigar(String cigar){
		this.cigar = cigar;
	}
	
	private void setMobileElementName(String name){
		this.mobile_element = name;
	}
	
	private void setNrOfMappings(SAMRecord rec, String tool) throws UnknownParamException{
		NrMappingsSAMRecordHolder recordHolder = NrMappingsSAMRecordHolderFactory.makeNrMappingsSAMRecordHolder(rec, tool);
		this.number_hits = Integer.toString(recordHolder.getNumberOfMappings());
	}
	
	private void setNrOfMappings(int nrMappings){
		this.number_hits = Integer.toString(nrMappings);
	}
	
	public void setHomoPolymer(String polymer){
		this.homopolymer = polymer;
	}
	
}
