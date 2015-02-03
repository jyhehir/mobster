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
	
	public void parse(String MEValueAttribute){
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
	
	public void parseMobileCategories(String[] categories){
		//do not do any checking anymore on the validness of the category when calling this function
		for (String category : categories){
			this.mobile_categories.add(category);
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


	public boolean addMobileCategoryByReference(String reference) throws InvalidCategoryException{
		
		String[] referenceSplit;
		String mobileCategory;
		
		//First check if the reference sequence contains substring OTHER_MOBILE_CATEGORY_ATTRIBUTE;
		//This allows for ME categorization outside the default Alu, L1, SVA, HERV categories.
		if (reference.contains(MobileDefinitions.OTHER_MOBILE_CATEGORY_ATTRIBUTE)){
			referenceSplit = reference.split(MobileDefinitions.OTHER_MOBILE_CATEGORY_ATTRIBUTE, -1);
			mobileCategory = referenceSplit[referenceSplit.length - 1].trim();
			if ( mobileCategory.length() == 0 ) {
				throw new InvalidCategoryException (reference + "can not be categorized\n." +  
			 " Tried to split on " + MobileDefinitions.OTHER_MOBILE_CATEGORY_ATTRIBUTE);
			}
			if (this.mobile_categories.contains(mobileCategory)){
				return false;
			}else{
				return this.mobile_categories.add(mobileCategory);	
			}

		}
		//If reference does not contain substring OTHER_MOBILE_CATEGORY_ATTRIBUTE
		else{
			for (String validCategory : MobileDefinitions.MOBILE_CATEGORIES){
				if (this.mobile_categories.contains(validCategory)){
					return false;
				}else if(reference.toUpperCase().startsWith(validCategory.toUpperCase())){
					return this.mobile_categories.add(validCategory);
				}
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
