package org.umcn.me.util;

import java.io.File;

import org.umcn.me.sam.InvalidCategoryException;
import org.umcn.me.sam.MobileSAMTag;


import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;


public class SAMReadNameClipper {
	
	public static void main (String[] args){
		File inFile = new File("D:/projects/mobile_finder_test/test_set/prealpha_splitread_pairedend/splitrefac_9mmNOsplit_anchors_csorted.bam");
		//clipSAMAndWrite(inFile, new File("D:/projects/mobile_finder_test/test_set/activemob100excl_40cov_bwatestRNcorr.bam"));
		getReadDetails(inFile);
	}
	
	public static void getReadDetails(File inFile){

		SAMFileReader samReader = new SAMFileReader(inFile);
		samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		for (SAMRecord rec : samReader){
			MobileSAMTag tag = new MobileSAMTag();
			System.out.println(rec.getAttribute("ME").toString());
			try {
				tag.parse(rec.getAttribute("ME").toString());
				System.out.println(tag.getMobileCigar());
				System.out.println(tag.getMobileElementName());
				System.out.println(tag.getNrOfMappings());
				System.out.println(tag.getMobileCategoryNames());
			} catch (InvalidCategoryException e) {
				
			}
			
		}
		
		samReader.close();
		

	}
	

	
	public static void clipSAMAndWrite(File inFile, File outFile){
		
		SAMFileReader samReader = new SAMFileReader(inFile);
		samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReader.getFileHeader(),
	    		false, outFile);
		
		String readName;
		
		for (SAMRecord record : samReader){
			readName = record.getReadName();
			readName = readName.substring(0, readName.length() - 1);
			record.setReadName(readName);
			samWriter.addAlignment(record);
			
			
			
		}
		
		samReader.close();
		samWriter.close();
		
	}

}
