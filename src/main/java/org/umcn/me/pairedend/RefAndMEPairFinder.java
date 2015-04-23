package org.umcn.me.pairedend;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.Vector;
import java.util.Map;
import java.util.HashMap;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.apache.commons.cli.*;
import org.umcn.gen.sequence.InvalidSequenceException;
import org.umcn.gen.sequence.Sequence;
import org.umcn.me.sam.InvalidCategoryException;
import org.umcn.me.sam.MobileSAMTag;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.samexternal.SAMWriting;
import org.umcn.me.samexternal.UnknownParamException;
import org.umcn.me.util.MobileDefinitions;

/*
 * IMPORTANT NOTE: STILL HAVE TO CHECK IF UPDATEREADSMAPPINGTOMULTIPLEME
 * WORKS PROPERLY AS BEST CASE BAM HAS NO READS THAT ALIGN AMBIGUOUSLY
 * ALSO A105A GONL SAMPLE HAD NO READS MAPPING AMBIGUOUSLY TO MOBILE REF
 */
/**
 * This class extracts all anchors for mates that map to the mobile reference
 * and meet the following conditions:
 * 
 * - For unique-unique read pairs an anchor is only valid and will be extracted
 * if one end maps to the mobile ref and the other does not. If both reads map
 * to the mobile ref or if both ends do not map to mobile ref, anchors will
 * not be extracted.
 * - For unique-multiple, multiple-unique and orphan-unique read pairs an anchor will be
 * extracted if the multiple end or orphan end maps at least once to the mobile ref.
 * 
 * For all anchors information about the mobile mate is added to the
 * ME attribute in the outputted BAM File. This attribute contains the following info:
 * - name of mobile element with best alignment
 * - cigar string of best alignment with mobile element
 * - The names of mobile elements superfamilies the read maps to. This value can be
 * used to detect whether read maps ambiguosly to multiple mobile superfamilies.
 * - Number of mobile elements the read maps to.
 * 
 * At the end all anchors are written to output bam file (*_anchors.bam).
 *
 * @update 3-7-2012: added HERV to mobile_elements categories
 * 
 * @author Djie
 */
public class RefAndMEPairFinder {

	public static Logger logger = Logger.getLogger("RefAndMEPairFinder"); 
	private static int MIN_POLYA_LEN = 9;
	private static int MAX_POLYA_MM = 1;
	private static int MAX_MAPPINGS = -1; // -1 do not use max mobiome mappings
	
	private static int MEMORY;
	private static String TMP;
	
	public static void main(String[] args) {
		Options options;
		File single = null;
		File multiple = null;
		File filtered = null;
		String out;
		String sampleName = "";
		String tool = "";
		Boolean paired;
		Properties prop = new Properties();
		
		
		HelpFormatter formatter = new HelpFormatter();
		BasicConfigurator.configure();
		
		options = createCmdOptions();
		
		if(args.length == 0){
			formatter.printHelp("java -Xmx4g -jar RefAndMEPairFinder.jar" , options);
		}else{
			CommandLineParser parser = new GnuParser();
			try {
				long start = System.currentTimeMillis();
				
				CommandLine line = parser.parse(options, args);
				
				if(line.hasOption("properties")){
					prop.load(new FileInputStream(new File(line.getOptionValue("properties"))));
					single = new File(prop.getProperty(MobileDefinitions.INFILE_FROM_MOBIOME_MAPPING).trim());
					//we are skipping the multiple input
					
					filtered = new File(prop.getProperty(MobileDefinitions.INFILE_FROM_POTENTIAL_MEI_FINDER).trim());
					out = prop.getProperty(MobileDefinitions.OUTFILE).trim();
					sampleName = prop.getProperty(MobileDefinitions.SAMPLE_NAME);
					tool = prop.getProperty(MobileDefinitions.MOBIOME_MAPPING_TOOL);
					paired = Boolean.parseBoolean(prop.getProperty(MobileDefinitions.PAIRED_END).trim());
					MIN_POLYA_LEN = Integer.parseInt(prop.getProperty(MobileDefinitions.POLY_A_LENGTH).trim());
					MAX_POLYA_MM = Integer.parseInt(prop.getProperty(MobileDefinitions.POLY_A_MAX_MISMATCHES).trim());
					
					if (prop.containsKey(MobileDefinitions.TMP)){
						TMP = prop.getProperty(MobileDefinitions.TMP).trim();
					}else{
						TMP = System.getProperty("java.io.tmpdir");
					}
					if (prop.containsKey(MobileDefinitions.GRIPS_MAX_REFSEQ_MAPPING)){
						MAX_MAPPINGS = Integer.parseInt(prop.getProperty(MobileDefinitions.GRIPS_MAX_REFSEQ_MAPPING));
					}
					
					MEMORY = Integer.parseInt(prop.getProperty(MobileDefinitions.MEMORY).trim());
					
					
				}else{
					single = new File(line.getOptionValue("single"));
					if(line.hasOption("multiple")){
						multiple = new File(line.getOptionValue("multiple"));
						logger.info("Checking for reads mapping to multiple superfamilies");
					}else{
						logger.info("Not checking for reads mapping to multiple superfamilies");
					}
					filtered = new File(line.getOptionValue("potential"));
					out = line.getOptionValue("out");
					tool = line.getOptionValue("tool", SAMDefinitions.MAPPING_TOOL_MOSAIK);
					sampleName = line.getOptionValue("samplename", sampleName);
					paired = line.hasOption("p");
					MAX_MAPPINGS = Integer.parseInt(line.getOptionValue("max_mapping", Integer.toString(MAX_MAPPINGS)));
					TMP = line.getOptionValue("tmp", System.getProperty("java.io.tmpdir"));
					MEMORY = Integer.parseInt(line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));
				}
				logger.info("Using potential .bam file: " + filtered);
				logger.info("Using output prefix: " + out);
				logger.info("Mapping tool used for mobile mapping: " + tool);
				logger.info("Paired-end reads: " + paired);
				logger.info("Setting sample name to: " + sampleName);
				
				runRefAndMePairFinder(single, multiple, filtered, out, tool, paired, sampleName);
				
				long end = System.currentTimeMillis();
				long millis = end - start;
				String time = String.format("%d min, %d sec", 
					    TimeUnit.MILLISECONDS.toMinutes(millis),
					    TimeUnit.MILLISECONDS.toSeconds(millis) - 
					    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
					);
				
				logger.info("RefAndMEPairFinder ran in : " + time);
				
			} catch (ParseException e) {
				logger.error("Error in parsing CLI arguments: " + e.getMessage());
			} catch (IOException ex){
				logger.error("Error in finding files: " + ex.getMessage());
			}
		}
 
	}
	
	//TODO remove the duplicated code associated with this code
	public static void runFromPropertiesFile(Properties prop) throws IOException{
		
		File single = null;
		File multiple = null;
		File filtered = null;
		String out;
		String sampleName = "";
		String tool = "";
		Boolean paired;
		
		long start = System.currentTimeMillis();
		
		single = new File(prop.getProperty(MobileDefinitions.INFILE_FROM_MOBIOME_MAPPING).trim());
		//we are skipping the multiple input
		
		filtered = new File(prop.getProperty(MobileDefinitions.INFILE_FROM_POTENTIAL_MEI_FINDER).trim());
		out = prop.getProperty(MobileDefinitions.OUTFILE).trim();
		sampleName = prop.getProperty(MobileDefinitions.SAMPLE_NAME);
		tool = prop.getProperty(MobileDefinitions.MOBIOME_MAPPING_TOOL);
		paired = Boolean.parseBoolean(prop.getProperty(MobileDefinitions.PAIRED_END).trim());
		MIN_POLYA_LEN = Integer.parseInt(prop.getProperty(MobileDefinitions.POLY_A_LENGTH).trim());
		MAX_POLYA_MM = Integer.parseInt(prop.getProperty(MobileDefinitions.POLY_A_MAX_MISMATCHES).trim());
		
		if (prop.containsKey(MobileDefinitions.GRIPS_MAX_REFSEQ_MAPPING)){
			MAX_MAPPINGS = Integer.parseInt(prop.getProperty(MobileDefinitions.GRIPS_MAX_REFSEQ_MAPPING));
		}
		
		if (prop.containsKey(MobileDefinitions.TMP)){
			TMP = prop.getProperty(MobileDefinitions.TMP).trim() + File.separator + "mob_" + Long.toString(System.nanoTime());			
		}else{
			TMP = System.getProperty("java.io.tmpdir") + File.separator + "mob_" + Long.toString(System.nanoTime());
		}
		
		File tmp = new File(TMP);
		
		if ( ! tmp.mkdir() ){
			throw new IOException("Can not create tmp directory: " + tmp);
		}
		MEMORY = Integer.parseInt(prop.getProperty(MobileDefinitions.MEMORY).trim());
		
		logger.info("Using potential .bam file: " + filtered);
		logger.info("Using output prefix: " + out);
		logger.info("Mapping tool used for mobile mapping: " + tool);
		logger.info("Paired-end reads: " + paired);
		logger.info("Setting sample name to: " + sampleName);
		logger.info("Using temp: " + tmp);
		
		try {
			runRefAndMePairFinder(single, multiple, filtered, out, tool, paired, sampleName);
			long end = System.currentTimeMillis();
			long millis = end - start;
			String time = String.format("%d min, %d sec", 
				    TimeUnit.MILLISECONDS.toMinutes(millis),
				    TimeUnit.MILLISECONDS.toSeconds(millis) - 
				    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
				);
			
			logger.info("RefAndMEPairFinder ran in : " + time);
		} catch (IOException e) {
			logger.error("RefAndMEPairFinder -> Error in finding files: " + e.getMessage());
		} finally{
			if (tmp != null && ! tmp.delete() ){
				logger.error("RefAndMEPairFinder -> Could not delete temp: " + tmp);
			}
		}
	}
	
	/**
	 * The do work function
	 * @param singleBam: bam file containing max 1 alignment per read against mobile ref.
	 * Should be readname sorted if paired data.
	 * @param multipleBam: bam file containing all alignments for all reads mapping
	 * multiple times against mobile ref. Should be readname sorted. (MOSAIK can output such a bam file)
	 * @param oriBam: bam file before mapping against the mobile reference, but after filtering (containing potential pairs)
	 * @param output: output bam file containing only anchors for mobile elements.
	 * @param tool: tool used for mobile reference mapping
	 * @param paired: true if paired-end reads are used
	 * @throws IOException
	 */
	public static void runRefAndMePairFinder(File singleBam, File multipleBam, File oriBam,
			 String output, String tool, Boolean paired, String sampleName) throws IOException{
			
		//TODO modify runRefAndMePairFinder to accept tool as a parameter, so different
		//mapping tools can be used more easily with this class. Now it can still be done
		//but multiple methods from this class need to be called seperately.
		Vector<String> exclusionReads = new Vector<String>();
		Map<String, MobileSAMTag> meReads;
		
		File nameSortedSingleBam = new File(singleBam.toString().replaceAll(".bam$",""));
		
		SAMWriting.writeSortedSAMorBAM(singleBam, nameSortedSingleBam, new File(TMP), MEMORY, SortOrder.queryname);
		//SAMWriting.writeNameSortedSAMorBAM(singleBam, nameSortedSingleBam, new File(TMP), MEMORY);
		
		exclusionReads = getUUReadsMappingBothToME(nameSortedSingleBam);
		
		meReads = getReadsMappingToME(nameSortedSingleBam, exclusionReads, tool);
		
		if(multipleBam != null){
			updateReadsMappingToMultipleME(multipleBam, meReads, exclusionReads);
		}
		
		writeRefCoordinateFile(oriBam, meReads, output, paired, sampleName);
		
	}

	/**
	 * This function returns all read names of unique-unique pairs where
	 * both reads map at least once to mobile reference
	 * NOTE: when there is a read name in the bam file called "dummy",
	 * this function may partially fail because of dumb coding ;-) see code
	 * and particular the initial value of var previousReadName
	 * @param bam name sorted bam file with mappings against mobile reference
	 * @return Vector of readnames of UU pairs with both mapping to mobile element
	 */
	public static Vector<String> getUUReadsMappingBothToME(File bam){
		SAMFileReader inputSam = new SAMSilentReader(bam);
		String readName = "";
		String previousReadName = "dummy";
		SAMRecord previousReadSam = null;
		Vector<String> exclReads = new Vector<String>();
		
		int c = 0;
		
		 for (SAMRecord samRecord : inputSam){
			 readName = samRecord.getReadName();
			 
			 //readName.length() - 1 to get rid of the last character in the read name, which with paired-end reads
			 //always end in 1 or 2.
			 if (readName.substring(0, readName.length() - 1).equals(previousReadName.substring(0, previousReadName.length() - 1))){
				 if(!samRecord.getReadUnmappedFlag() && !previousReadSam.getReadUnmappedFlag()){
					 exclReads.add(readName);
					 exclReads.add(previousReadName);
					 c++;
				 }
	
			 }
			 previousReadName = samRecord.getReadName();
			 previousReadSam = samRecord;
		 }
		
		logger.info("Found: " + c + " unique-unique pairs with both ends mapping to mobile ref, excluding filtered reads");
		logger.info("These reads are filtered");
		 
		inputSam.close();
		
		return exclReads;
	}

	/**
	 * This function puts all reads which map at least once to a mobile element
	 * inside a HashMap alongside additional information (to which mobile element
	 * the read maps, cigar string and number of hits to mobile elements).
	 * The read may not exist in the exclusion vector, otherwise the read will
	 * not be included in the hashmap
	 * @param bam file with reads mapped to mobile reference (NAME SORTED!)
	 * @param exclusion vector containing read names which will not be put in
	 * the hashmap, even when the read maps to a mobile element.
	 * @param name of tool used for mapping against mobile reference (bwa for bwa mapping)
	 * @return HashMap containing reads and additional information about read
	 * and mapping info of reads mapping at least once to mobile element and
	 * not occuring in exclusion vector.
	 */
	public static Map<String, MobileSAMTag> getReadsMappingToME(File bam,
			Vector<String> exclusion, String tool){
		
		Map<String, MobileSAMTag> meReads = new HashMap<String, MobileSAMTag>();
		SAMFileReader inputSam = new SAMSilentReader(bam);	
		int mobileReadCounter = 0;
		int skippedReadsBecauseOfMultiMappings = 0;
		
		logger.info("Using max mobiome mappings of (-1 is disabled) : " + MAX_MAPPINGS);
		
		for (SAMRecord samRecord : inputSam){
			
			String recordName = samRecord.getReadName();
			boolean splitRead = recordName.startsWith(SAMDefinitions.SPLIT_MAPPING);
			
			if(MAX_MAPPINGS != -1){
				MobileSAMTag tempTag = new MobileSAMTag();
				try {
					tempTag.build(samRecord, tool);
					if (tempTag.getNrOfMappings() > MAX_MAPPINGS){
						skippedReadsBecauseOfMultiMappings++;
						continue;
					}
				} catch (InvalidCategoryException e) {
					logger.error(e.getMessage());
				} catch (UnknownParamException e) {
					logger.error(e.getMessage());
				}
				
			}
			
			if (!samRecord.getReadUnmappedFlag() && !exclusion.contains(samRecord.getReadName())){
				mobileReadCounter++;
				
				MobileSAMTag mobileTag = new MobileSAMTag();
				try {
					if(splitRead){
						Sequence seq = new Sequence(samRecord.getReadString());
						if(samRecord.getReadNegativeStrandFlag()){
							seq.reverseComplement();
						}
						String homoPolymer = "";
						if (recordName.charAt(1) == SAMDefinitions.LEFT_CLIPPED){
							homoPolymer = seq.getPolyAOrTMapping(MIN_POLYA_LEN, MAX_POLYA_MM, false);
						}else if(recordName.charAt(1) == SAMDefinitions.RIGHT_CLIPPED){
							homoPolymer = seq.getPolyAOrTMapping(MIN_POLYA_LEN, MAX_POLYA_MM, true);
						}
						if("polyA".equalsIgnoreCase(homoPolymer) || "polyT".equalsIgnoreCase(homoPolymer)){
							mobileTag.setHomoPolymer(homoPolymer);
						}
					}
					mobileTag.build(samRecord, tool);
					meReads.put(samRecord.getReadName(), mobileTag);
				} catch (InvalidCategoryException e) {
					logger.error(e.getMessage());
				} catch (UnknownParamException e) {
					logger.error(e.getMessage());
				} catch (InvalidSequenceException e){
					logger.error(e.getMessage());
				}
			}else if(splitRead && samRecord.getReadUnmappedFlag()){
				Sequence splitSeq;
				try {
					String homoPolymer = "";
					splitSeq = new Sequence(samRecord.getReadString());
					if (recordName.charAt(1) == SAMDefinitions.LEFT_CLIPPED){
						homoPolymer = splitSeq.getPolyAOrTMapping(MIN_POLYA_LEN, MAX_POLYA_MM, false);
					}else if(recordName.charAt(1) == SAMDefinitions.RIGHT_CLIPPED){
						homoPolymer = splitSeq.getPolyAOrTMapping(MIN_POLYA_LEN, MAX_POLYA_MM, true);
					}
					if("polyA".equalsIgnoreCase(homoPolymer) || "polyT".equalsIgnoreCase(homoPolymer)){
						MobileSAMTag mobileTag = new MobileSAMTag();
						mobileTag.setHomoPolymer(homoPolymer);
						mobileTag.build(samRecord, tool);
						meReads.put(samRecord.getReadName(), mobileTag);
					}
					
				} catch (InvalidSequenceException e) {
					logger.error(e.getMessage());
				} catch (InvalidCategoryException e){
					logger.error(e.getMessage());
				}catch (UnknownParamException e) {
					logger.error(e.getMessage());
				}
			}
		}
		logger.info(mobileReadCounter + " reads map to at least 1 mobile element");
		logger.info("Skipped reads because of mobiome multimappings: " + skippedReadsBecauseOfMultiMappings);
		inputSam.close();
		return meReads;
	}

	/**
	 * NOTE THIS FUNCTION IS STILL PARTLY UNTESTED!
	 * In theory this function should update the HashMap containing reads which map to mobile
	 * reference. Reads which map to more than one mobile element superfamily will get
	 * information regarding this updated inside the hashmap (mecats slot).
	 * @param bam multiple bam file, name sorted from MOSAIK
	 * @param meReads HashMap containing read names as keys and MobileSAMTag as values
	 * containing mobile read mapping info
	 * @param exclusion vector of reads not to include in this analysis.
	 */
	public static void updateReadsMappingToMultipleME(File bam, Map<String, MobileSAMTag> meReads,
			Vector<String> exclusion){
		
		SAMFileReader inputSam = new SAMSilentReader(bam);
		
		String readName;
		int c = 0;
				
		for (SAMRecord samRecord : inputSam){
			readName = samRecord.getReadName();
			if(!exclusion.contains(readName) && meReads.containsKey(readName)){
				try {
					if (meReads.get(readName).getHomoPolymer().equals("") &&
							meReads.get(readName).addMobileCategoryByReference(samRecord.getReferenceName())){
						logger.info(readName);
						c++;
					}
				} catch (InvalidCategoryException e) {
					logger.error(e.getMessage());
				}
			}
		}
		logger.info(c + " reads align to more than 1 Mobile Element Family");
		inputSam.close();
		
	}

	/**
	 * This function writes the bam file with only reads functioning as
	 * anchor points for MEI events.
	 * @param originalBam: bam file with reads after filtering but before
	 * mapping to mobile reference (bam file containing potential pairs).
	 * @param meReads: hashmap containing reads mapping to mobile elements and additional info as a MobileSAMTag.
	 * Read Name acts as key.
	 * @param outputBam: name for anchor only output bam file
	 * @param paired: whether reads are paired-end (true)
	 * @throws IOException
	 */
	public static void writeRefCoordinateFile(File originalBam,
			Map<String, MobileSAMTag> meReads, String out, Boolean paired, String sampleName) throws IOException{
		
		String mobileSAMValue = "";
		StringBuilder anchor = new StringBuilder();
		boolean splitRead;
		
		SAMFileReader inputSam = new SAMSilentReader(originalBam);
		
		SAMFileWriter mateClusterSam = SAMWriting.makeSAMWriter(new File(out + "_anchors.bam"),
				inputSam.getFileHeader(), new File(TMP), MEMORY, SAMFileHeader.SortOrder.coordinate);
		
		SAMFileWriter splitClusterSam = SAMWriting.makeSAMWriter(new File(out + "_splitanchors.bam"),
				inputSam.getFileHeader(), new File(TMP), MEMORY, SAMFileHeader.SortOrder.coordinate);
		
		int c = 0;
		int d = 0;
		for (SAMRecord samRecord : inputSam){
			anchor.append(samRecord.getReadName());
			splitRead = false;
			if(paired){
				anchor.append(SAMDefinitions.READ_NUMBER_SEPERATOR);
				splitRead = anchor.toString().startsWith(SAMDefinitions.SPLIT_MAPPING);
				if(samRecord.getFirstOfPairFlag()){
					if(splitRead){
						anchor.append("1");
					}else{
						anchor.append("2");
					}
				}else{
					if(splitRead){
						anchor.append("2");
					}else{
						anchor.append("1");
					}
				}
			}
			if(meReads.containsKey(anchor.toString())){
				//in this case samRecord is anchor point for mobile Mate
				mobileSAMValue = meReads.get(anchor.toString()).toString();

				samRecord.setAttribute(MobileDefinitions.SAM_TAG_MOBILE, mobileSAMValue);
				samRecord.setAttribute(MobileDefinitions.SAM_TAG_SAMPLENAME, sampleName);

				
				if(splitRead){
					splitClusterSam.addAlignment(samRecord);
					c++;
				}else{
					mateClusterSam.addAlignment(samRecord);
					d++;
				}
				
			}
			anchor.setLength(0);
		}
		logger.info("Nr of split reads written to  split anchor file: " + c);
		logger.info("Nr of cluster mates written to anchor file: " + d);
		inputSam.close();
		mateClusterSam.close();
		splitClusterSam.close();
	}

	public static Options createCmdOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("BAM file (NAME-sorted!) from mobile ref mapping," +
				"containing at most 1 mapping per read");
		
		options.addOption(OptionBuilder.create("single"));
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("BAM file (NAME-sorted!) from mobile ref mapping,"
				+ "containing all mappings of multiply mapped reads");
		
		options.addOption(OptionBuilder.create("multiple"));
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("BAM file containing pairs potentially supporting MEI events" +
				" (output BAM file of PotentialMEIFinder");
		
		options.addOption(OptionBuilder.create("potential"));
		
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output prefix for BAM file containing anchors");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("Mapping tool");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Mapping tool used for mobile mapping, default: " + 
				SAMDefinitions.MAPPING_TOOL_MOSAIK);
		
		options.addOption(OptionBuilder.create("tool"));
		
		OptionBuilder.withArgName("Paired-End Data");
		OptionBuilder.withDescription("Whether the original mapping is done using paired-end reads." +
				" Leave -p out when using single-end fragment data.");
		
		options.addOption(OptionBuilder.create("p"));
		
		OptionBuilder.withArgName("SampleName");
		OptionBuilder.withDescription("name of sample");
		options.addOption(OptionBuilder.create("samplename"));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Tmp directory to use, when writing large files");
		
		options.addOption(OptionBuilder.create("tmp"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Alter the number of records in memory (increase to decrease number of file handles). Default: " + SAMWriting.MAX_RECORDS_IN_RAM);
		
		options.addOption(OptionBuilder.create("max_memory"));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Use property file instead of command line arguments, command line arguments will be skipt. Values in property file will be used");
		
		options.addOption(OptionBuilder.create("properties"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum number of mappings a read may map to the mobiome");
		options.addOption(OptionBuilder.create("max_mapping"));
		
		return options;
	}
	
}
