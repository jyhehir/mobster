package org.umcn.me.pairedend;

import net.sf.samtools.*;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import java.io.*;
import java.util.*;
import org.umcn.gen.parsers.FastaParserFast;
import org.umcn.gen.sam.SAMDefinitions;
import org.umcn.gen.sam.UnknownParamException;



/**
 * This class is used to find read pairs that could potentially indicate
 * a MEI event:
 * - Read pairs which map for one end uniquely to the genome and for
 *   the other end multiple times to the genome; but only when
 *   the reads are not properly paired**.
 * - Read pairs which map for for one end uniquely to the genome and
 *   for the other end 0 times to the genome (orphans).
 * - Read pairs which map for both ends uniquely to genome but
 *   are not properly paired**.
 *   
 *   **unproperly paired:
 *   	- both reads map on + strand
 *      - both reads map on - strand
 *      - read mapping on + strand has higher genomic coordinate
 *        than read mapping on - strand
 *      - proper FR orientation but insert size not in
 *        mean fragment size +/- local search
 *        
 *  All potential MEI supporting read pairs are written to a bam file.
 *  All reads in a read pair which may potentially map (orphans, multiple mappings,
 *  unproper pairs for unique-unique pairs) to mobile references
 *  are written to seperate fq files.
 * @author Djie
 *
 */
public class PotentialMEIFinder {
	
	public static Logger logger = Logger.getLogger("PotentialMEIFinder");
	
	private String pair_number_seperator;

	private String mapping_tool;

	private String output_prefix;

	private String bam_mapped_to_ref;

	
	
	public PotentialMEIFinder(String pairNumberSeperator, String inBam, String outPrefix,
			String tool) throws UnknownParamException{
		BasicConfigurator.configure();
		
		this.pair_number_seperator = pairNumberSeperator;
		this.bam_mapped_to_ref = inBam;
		this.output_prefix = outPrefix;
		this.mapping_tool = tool;
		
		supportedParamCheck(SAMDefinitions.MAPPING_TOOLS, tool);
		
	}
	
	/**
	 * Checks for each sam/bam entry whether a
	 * pair is potentially supporting a MEI event. Potential pairs are written to
	 * a BAM file with the prefix specified by the out prefix private field appended with _potential.bam.
	 * Furthermore potential mates are written to one .fq file (out_prefix + "_potential.fq")
	 * 
	 * @throws IOException 
	 * @throws UnknownParamException 
	 */
	public void extractPotentialMobileReads() throws IOException {
	
		long start = System.nanoTime();
		File input = new File(this.bam_mapped_to_ref);
	    int mappingsRead; //Nr mappings for current SAM Record Read
	    int mappingsMate; //Nr mappings for mate of current SAM Record Read
	    int nrPotentialMobileReads = 0;
	    int nrMatesOfPotentialMobileReads = 0;
	    int c = 0;
	
	    File outputBam = new File(this.output_prefix + "_potential.bam");
	    PrintWriter outFq = new PrintWriter(new FileWriter(this.output_prefix + "_potential.fq"), true);      
	
	    SAMFileReader inputSam = new SAMFileReader(input);
	    inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	    
	    SAMFileReader inputSam2 = new SAMFileReader(input);
	    inputSam2.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
	    
	    SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
	    		false, outputBam);
		
	    for (SAMRecord samRecord : inputSam) {
	    	c++;
	    	
	    	if (c % 10000 == 0){
	    		System.out.println(c);
	    	}
	    	
	    	mappingsRead = getNumberOfMappings(samRecord, this.mapping_tool);
	    	mappingsMate = getNumberOfMappingsOfMate(samRecord, this.mapping_tool, inputSam2);
			
			if (isReadPotentialMobile(mappingsRead, mappingsMate, samRecord)){
				FastaParserFast.writeFastQToStream(outFq, makeFastQHeaderFromSAMRecord(samRecord, this.pair_number_seperator),
						samRecord.getReadString(), samRecord.getBaseQualityString());
				outputSam.addAlignment(samRecord);
				nrPotentialMobileReads++;
			}else if (isReadMateOfPotentialMobile(mappingsRead, mappingsMate, samRecord)){
				outputSam.addAlignment(samRecord);
				nrMatesOfPotentialMobileReads++;
			}
	    }
	
	    logger.info("Nr of found potential mobile reads: " + nrPotentialMobileReads);
	    logger.info("Nr of found mates of potential reads: " + nrMatesOfPotentialMobileReads);
	
	    inputSam.close();
	    outputSam.close();
	    inputSam2.close();
	    outFq.close();
	    long end = System.nanoTime();
	    System.out.println((double) (end - start) / 1000000000.0);
	}

	private static boolean isReadPotentialMobile(int mappingsRead,
			int mappingsMate, SAMRecord samRecord) {
		
		boolean found = false;
		
		if(isFU(mappingsRead, mappingsMate)){
			found = true;
		}else if(isMU(samRecord, mappingsRead, mappingsMate)){
			found = true;
		}else if(isUnproperUU(samRecord, mappingsRead, mappingsMate)){
			found = true;
		}
		return found;
	}
	
	private static boolean isReadMateOfPotentialMobile(int mappingsRead,
			int mappingsMate, SAMRecord samRecord) {
		
		boolean found = false;
		
		if(isUF(mappingsRead, mappingsMate)){
			found = true;
		}else if(isUM(samRecord, mappingsRead, mappingsMate)){
			found = true;
		}else if(isUnproperUU(samRecord, mappingsRead, mappingsMate)){
			found = true;
		}
		return found;
		
	}
    
    private static String makeFastQHeaderFromSAMRecord(SAMRecord samRecord, String sep) {
		String pairNumber = "1";
		StringBuilder header = new StringBuilder();
		
		header.append("@");
		header.append(samRecord.getReadName());
		header.append(sep);
		
		if (!samRecord.getFirstOfPairFlag()) pairNumber = "2";
		
		header.append(pairNumber);
		
		return header.toString();
			
	}

	private int getNumberOfMappingsOfMate(SAMRecord rec,
			String tool, SAMFileReader samReader) {

    	int matches = 0;
    	
    	if (SAMDefinitions.MAPPING_TOOL_MOSAIK.equals(tool)){
            matches = getNumberOfMappingsMateMosaik(rec);
           
    	}else if (SAMDefinitions.MAPPING_TOOL_BWA.equals(tool)){
    		matches = getNumberOfMappingsMateBWA(rec, samReader);
    	}
    	
    	return matches;
	}

	private static int getNumberOfMappingsMateMosaik(SAMRecord rec) {
		int nrMappings;
		String pairMappingInfo; //ZA tag outputted by MOSAIK
		String singleReadMappingInfo;
		
 
		if (rec.getAttribute(SAMDefinitions.MOSAIK_MAPPINGINFO_ATTRIBUTE) != null) {
			pairMappingInfo = rec.getAttribute(SAMDefinitions.MOSAIK_MAPPINGINFO_ATTRIBUTE).toString();
			if (rec.getFirstOfPairFlag()) {
				singleReadMappingInfo = pairMappingInfo.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ2];
			} else {
				singleReadMappingInfo = pairMappingInfo.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ1];
			}
			nrMappings = Integer.parseInt(singleReadMappingInfo.split(";")[SAMDefinitions.MOSAIK_MAPPINGINFO_NRHITS]);
		}else{
			nrMappings = 0;
		}
		return nrMappings;
	}

	public static int getNumberOfMappings(SAMRecord rec, String tool){
    	
    	int matches = 0;
    	
    	if (SAMDefinitions.MAPPING_TOOL_MOSAIK.equals(tool)){
            matches = getNumberOfMappingsMosaik(rec);       
    	}else if (SAMDefinitions.MAPPING_TOOL_BWA.equals(tool)){
    		matches = getNumberOfMappingsBWA(rec);
    	}
    	
    	return matches;
    }

	private static int getNumberOfMappingsBWA(SAMRecord rec) {
		
		int nrMappings;
		
		if (rec.getReadUnmappedFlag()){
			nrMappings = 0;
		}else{
			if(rec.getAttribute(SAMDefinitions.BWA_NR_OPTIMAL_HITS_ATTRIBUTE) != null){
				nrMappings = Integer.parseInt(
						rec.getAttribute(SAMDefinitions.BWA_NR_OPTIMAL_HITS_ATTRIBUTE).toString());
			}else{
				nrMappings = 1;
			}

		}
		
		return nrMappings;
	}
	
	//TODO finish implementation
	private int getNumberOfMappingsMateBWA(SAMRecord rec, SAMFileReader inputSam) {
		
		int nrMappings;
		
		if (rec.getMateUnmappedFlag()){
			nrMappings = 0;
		}else{
			SAMRecord mate = inputSam.queryMate(rec);
			if(mate.getAttribute(SAMDefinitions.BWA_NR_OPTIMAL_HITS_ATTRIBUTE) != null){
				nrMappings = Integer.parseInt(
						mate.getAttribute(SAMDefinitions.BWA_NR_OPTIMAL_HITS_ATTRIBUTE).toString());
			}else{
				nrMappings = 1;
			}
		}
		
		return nrMappings;
	}

	private static int getNumberOfMappingsMosaik(SAMRecord rec) {
		int nrMappings;
		String pairMappingInfo; //ZA tag outputted by MOSAIK
		String singleReadMappingInfo;
		
		if (rec.getAttribute(SAMDefinitions.MOSAIK_MAPPINGINFO_ATTRIBUTE) != null) {
			pairMappingInfo = rec.getAttribute(SAMDefinitions.MOSAIK_MAPPINGINFO_ATTRIBUTE).toString();
			if (rec.getFirstOfPairFlag()) {
				singleReadMappingInfo = pairMappingInfo.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ1];
			} else {
				singleReadMappingInfo = pairMappingInfo.split("><")[SAMDefinitions.MOSAIK_MAPPINGINFO_READ2];
			}
			nrMappings = Integer.parseInt(singleReadMappingInfo.split(";")[SAMDefinitions.MOSAIK_MAPPINGINFO_NRHITS]);
		}else{
			nrMappings = 0;
		}
		return nrMappings;
	}
    
    /**
     * This function extracts reads from a fastq file (fqFile) when the read name
     * in the idSet equals the read name in the fqFile. All extracted reads are
     * written in fastq format into outFile.
     * NOTE: IT IS ASSUMED THAT ENTRIES IN THE FASTQ FILE ALWAYS CONTAIN 4 LINES!
     * 
     * @param idSet Set containing read ids (should contain fullname with prepended "@" to readname)
     * @param fqFile String containing location of fastq file from which reads should be extracted
     * @param outFile String containing location + file name to which extracted reads should be written
     * @param append Boolean append to outFile or not
     * @return the number of read names that overlap in idSet and fqFile.
     * @throws IOException
     */
	@Deprecated
    public static int createFqSubset(Set<String> idSet, String fqFile, String outFile, boolean append) throws IOException{
    	PrintWriter out;
    	BufferedReader br;
    	String line;
    	boolean hit = false;
    	int lineCounter = 0;
    	int found = 0;
    	
    	br = new BufferedReader(new FileReader(fqFile));
    	out = new PrintWriter(new FileWriter(outFile, append), true);
    	
    	while ((line = br.readLine()) != null){
    		if(line.startsWith("@") && idSet.contains(line.trim())){
    			found += 1;
    			hit = true;
    			out.println(line.trim());
    		}else if(hit == true){
    			lineCounter++;
    			out.println(line.trim());
    			if(lineCounter == 3){
    				lineCounter = 0;
    				hit = false;
    			}
    			
    		}
    	}
    	br.close();
    	out.close();
    	return found;
    }
    
    /**
     * Is unproper unique-unique pair?
     * @param record
     * @param matches1 number of times read maps to reference
     * @param matches2 number of times mate of read maps to reference
     * @return
     */
    public static boolean isUnproperUU(SAMRecord record, int matches1, int matches2){
    	if(matches1 == 1 && matches2 == 1 && !record.getProperPairFlag()){
    		return true;
    	}
    	return false;
    }
    
    /**
     * Is first read orphan (FILTERED) and second read uniquely mapped?
     * @param hits1 number of times read maps to reference
     * @param hits2 number of times mate of read maps to reference
     * @return
     */
    public static boolean isFU(int hits1, int hits2){
    	if(hits1 == 0 && hits2 == 1){
    		return true;
    	}
    	return false;
    }
    
    /**
     * Is first read mapped multiple times and second read uniquely?
     * @param record
     * @param matches1 number of times read maps to reference
     * @param matches2 number of times mate of read maps to reference
     * @return
     */
    public static boolean isMU(SAMRecord record, int matches1, int  matches2){
    	if(matches1 > 1 && matches2 == 1 && !record.getProperPairFlag()){
    		return true;
    	}
    	return false;
    }
    
    /**
     * Is first read mapped uniquely and is second read orphan?
     * @param hits1 number of times read maps to reference
     * @param hits2 number of times mate of read maps to reference
     * @return
     */
    public static boolean isUF(int hits1, int hits2){
    	if(hits1 == 1 && hits2 == 0){
    		return true;
    	}
    	return false;
    }
    
    /**
     * Maps first read uniquely and second read multiple times?
     * @param record
     * @param matches1 number of times read maps to reference
     * @param matches2 number of times mate of read maps to reference
     * @return
     */
    public static boolean isUM(SAMRecord record, int matches1, int  matches2){
    	if(matches1 == 1 && matches2 > 1 && !record.getProperPairFlag()){
    		return true;
    	}
    	return false;
    }

	private static void supportedParamCheck(String[] params, String param) throws UnknownParamException{
		
		for (String p : params){
			if (p.equals(param)){
				return;
			}
		}
		
		throw new UnknownParamException("Unsupported argument: " + param);
	}
    
}
