package org.umcn.me.samexternal;

import java.io.File;
import java.util.Vector;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

public class MAPQFilterExplanation {

	private static String test_sam = "./test/test_sam.sam";
	
	public static void main(String[] args) {
		
		File samFile = new File(test_sam);
		int threshold = 20;
		
		mapqFilteringWayOne(samFile, threshold);
		
		mapqFilteringWayTwo(samFile);
		
		filteringUsingMultipleFilters(samFile);
	}
	
	//NOTE: preferred way when only filtering for one criterium
	public static void mapqFilteringWayOne(File sam, int mapqthreshold){
		SAMFileReader reader = new SAMFileReader(sam); //can also be a .bam file
		Vector<SAMRecord> approvedReads = new Vector<SAMRecord>();
		
		//loop through all reads in bam file
		for (SAMRecord record : reader){
			if (record.getMappingQuality() >= mapqthreshold){
				System.out.println("MAPQ above or equal: " + mapqthreshold);
				System.out.println(record.getSAMString());
				approvedReads.add(record);
			}else{
				System.out.println("MAPQ under: " + mapqthreshold);
				System.out.println(record.getSAMString());
			}
		}
		System.out.println("Number of approved reads way one: " + approvedReads.size());
		reader.close();
	}
	
	private static void mapqFilteringWayTwo(File sam) {
		
		SAMFileReader reader = new SAMFileReader(sam);
		//note BAMReadFilter can sometimes still be a bit quirky:
		//for instance when reads are also mapped against contig.
		//also filtering upon number of mappings per read will 
		//not work when you have not specified the used mapping tool
		BAMReadFilter filter = new BAMReadFilter();
		
		
		File outFile = new File("D:/out_test.bam");
		
		SAMFileHeader samFileHeader = reader.getFileHeader();
		System.out.println(samFileHeader.getTextHeader());
		
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader,
	    		false, outFile);
		//SAMFileWriter writer = new SAMFileWriter(outFile);
		filter.setMAPQFilter(20);
		Vector<SAMRecord> approvedReads = new Vector<SAMRecord>();
		
		for (SAMRecord record : reader){
			if(filter.passesFilter(record, false)){
				approvedReads.add(record);
				writer.addAlignment(record);
			}
		}
		
		writer.close();
		
		System.out.println("Number of approved reads way two: " + approvedReads.size());
		reader.close();
		
	}


	private static void filteringUsingMultipleFilters(File sam) {
		SAMFileReader reader = new SAMFileReader(sam);
		BAMReadFilter filter = new BAMReadFilter();
		Vector<SAMRecord> approvedReads = new Vector<SAMRecord>();
		filter.setAverageQualFilter(20.0);
		filter.setMaxEditDistance(1);
		
		for (SAMRecord record : reader){
			if(filter.passesFilter(record, false)){
				approvedReads.add(record);
			}
		}
		
		System.out.println("Number of approved reads using different filters: " + approvedReads.size());
		reader.close();
		
		
	}
}
