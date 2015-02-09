package org.umcn.me.pairedend;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.me.sam.PotentialMobilePairIterator;
import org.umcn.me.samexternal.NrMappingsSAMRecordHolder;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.samexternal.SAMRecordHolderPair;
import org.umcn.me.samexternal.SAMWriting;
import org.umcn.me.util.BAMCollection;
import org.umcn.me.util.MobileDefinitions;

import com.google.code.jyield.YieldUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;

public class PotentialMEIReadFinder {
	
	public static Logger logger = Logger.getLogger("PotentialMEIReadFinder");
	private static int min_avg_qual = 20;
	private static int min_anchor_mapq = 20;
	
	public static void main(String[] args) {
		
		BasicConfigurator.configure();
		
		Options options = addCmdOptions();
		HelpFormatter formatter = new HelpFormatter();	
		Properties props = new Properties();
		
		String infile;
		String outfile;
		String tool;
		String tmp;
		Boolean useSplit;
		int minClipping;
		int maxClipping;
		int memory;
		File propertiesFile;
		if(args.length == 0){
			formatter.printHelp("java -Xmx4g -jar PotentialMEIFinder.jar", options);
		}else{
			CommandLineParser parser = new GnuParser();
			try {
				long start = System.currentTimeMillis();
				
				CommandLine line = parser.parse(options, args);
				
				if (line.hasOption("properties")){
					propertiesFile = new File(line.getOptionValue("properties"));
					props.load(new FileInputStream(propertiesFile));
					infile = props.getProperty(MobileDefinitions.INFILE).trim();
					outfile = props.getProperty(MobileDefinitions.OUTFILE).trim();
					
					if (props.containsKey(MobileDefinitions.MAPPING_TOOL)){
						tool = props.getProperty(MobileDefinitions.MAPPING_TOOL).trim();
					}else{
						//3-12-2014:Changed from mosaik to unspecified
						tool = SAMDefinitions.MAPPING_TOOL_UNSPECIFIED;
					}
					useSplit = Boolean.parseBoolean(props.getProperty(MobileDefinitions.USE_SPLIT).trim());
					
					if (props.containsKey(MobileDefinitions.TMP)){
						tmp = props.getProperty(MobileDefinitions.TMP).trim();
					}else{
						tmp = System.getProperty("java.io.tmpdir");
					}
					
					minClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_CLIPPING).trim());
					maxClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_CLIPPING).trim());
					
					memory = Integer.parseInt(props.getProperty(MobileDefinitions.MEMORY).trim());
					
					min_avg_qual = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_QUAL).trim());
					if (props.containsKey(MobileDefinitions.MIN_MAPQ_ANCHOR)){
						min_anchor_mapq = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_MAPQ_ANCHOR).trim());
					}
					
					
				}else{
					infile = line.getOptionValue("in");
					outfile = line.getOptionValue("out");
					tool = line.getOptionValue("tool", SAMDefinitions.MAPPING_TOOL_UNSPECIFIED); //3-12-2014:Changed from mosaik to unspecified
					useSplit = line.hasOption("split");
					tmp = line.getOptionValue("tmp", System.getProperty("java.io.tmpdir"));
					minClipping = Integer.parseInt(line.getOptionValue("minclip", "35")); //30-1-2015: typo in option parsing min --> minclip
					maxClipping = Integer.parseInt(line.getOptionValue("maxclip", "7")); //30-1-2015: typo in option parsing max --> maxclip
					memory = Integer.parseInt(line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));
					min_avg_qual = Integer.parseInt(line.getOptionValue("min_avg_qual", Integer.toString(min_avg_qual)));
					min_anchor_mapq = Integer.parseInt(line.getOptionValue("mapq",Integer.toString(min_anchor_mapq)));
				}
				
				logger.info("Using infile: " + infile);
				logger.info("Using out prefix: " + outfile);
				logger.info("Used mapping tool: " + tool);
				logger.info("Including split reads: " + useSplit);
				logger.info("Maximum records in memory: " + memory);
				logger.info("Temp directory: " + tmp);
				logger.info("Split reads should have min avg qual of: " + min_avg_qual);
				
				if(useSplit){
					logger.info("Minimum clipping of: " + minClipping);
					logger.info("Maximum clipping of other side: " + maxClipping);
				}
				
				if (tool.equals(SAMDefinitions.MAPPING_TOOL_UNSPECIFIED)){
					logger.info("Using mapq of: " + min_anchor_mapq +  "for defining anchors");
				}
				
				//TODO1: Open up the .fq and .bam writer here, then run the potential MEIFinder
				
				//Loop over runPotentialMEIFinder depending on the number of BAM files given
				//runPotentialMEIFinder(infile, outfile, tool, tmp, useSplit, minClipping, maxClipping, memory); 
				
				
				//TODO3: Close the .fq and .bam writer here.
				
				long end = System.currentTimeMillis();
				long millis = end - start;
				String time = String.format("%d min, %d sec", 
					    TimeUnit.MILLISECONDS.toMinutes(millis),
					    TimeUnit.MILLISECONDS.toSeconds(millis) - 
					    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
					);
				
				logger.info("[PMRF] PotentialMEIReadFinder ran in : " + time);
			} catch (ParseException e) {
				logger.error("[PMRF] Error in parsing CLI arguments: " + e.getMessage());
			}  catch (FileNotFoundException e){
				logger.error("[PMRF] Could not create or find file / directory");
				logger.error(e.getMessage());
			} catch (IOException e) {
				logger.error("[PMRF] IO error");
				logger.error(e.getMessage());
			}
		}
		   
	}
	
	//TODO remove duplicate code for this function
	public static void runFromProperties(Properties props){
		
		//TODO
//		String infile;
		String outfile;
		String tool;
		String tmp;
		Boolean useSplit;
		int minClipping;
		int maxClipping;
		int memory;
		
		String[] samples;
		String[] bams;
		
		long start = System.currentTimeMillis();
		
		//TODO
		//infile = props.getProperty(MobileDefinitions.INFILE).trim();
		outfile = props.getProperty(MobileDefinitions.OUTFILE).trim();
		
		if (props.containsKey(MobileDefinitions.MAPPING_TOOL)){
			tool = props.getProperty(MobileDefinitions.MAPPING_TOOL).trim();
		}else{
			//3-12-2014:Changed from mosaik to unspecified
			tool = SAMDefinitions.MAPPING_TOOL_UNSPECIFIED;
		}
		
		useSplit = Boolean.parseBoolean(props.getProperty(MobileDefinitions.USE_SPLIT).trim());
		
		if (props.containsKey(MobileDefinitions.TMP)){
			tmp = props.getProperty(MobileDefinitions.TMP).trim();
		}else{
			tmp = System.getProperty("java.io.tmpdir");
		}
		
		//TODO: Sample parsing
		samples = props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0);
		bams = props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0);
		
		//Check to see whether all values are actually unique
		if (samples.length != new HashSet<String>(Arrays.asList(samples)).size() ||
				bams.length != new HashSet<String>(Arrays.asList(bams)).size()){
			logger.error("Supplied bams and/or supplied sample names are not unique");
			System.exit(1);
		}
		
		//TODO: If multiple samples were found, turn on the multiple_sample_calling option in the properties from Mobster.java
		
		minClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_CLIPPING).trim());
		maxClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_CLIPPING).trim());
		
		memory = Integer.parseInt(props.getProperty(MobileDefinitions.MEMORY).trim());
		
		min_avg_qual = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_QUAL).trim());
		
		logger.info("Using infile(s): " + props.getProperty(MobileDefinitions.INFILE));
		logger.info("Using out prefix: " + outfile);
		logger.info("Used mapping tool: " + tool);
		logger.info("Including split reads: " + useSplit);
		logger.info("Maximum records in memory: " + memory);
		logger.info("Temp directory: " + tmp);
		logger.info("Split reads should have min avg qual of: " + min_avg_qual);
		logger.info("Number of samples detected: " + samples.length);
		
		if(useSplit){
			logger.info("Minimum clipping of: " + minClipping);
			logger.info("Maximum clipping of other side: " + maxClipping);
		}
		
		if (tool.equals(SAMDefinitions.MAPPING_TOOL_UNSPECIFIED)){
			min_anchor_mapq = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_MAPQ_ANCHOR));
			logger.info("Using mapq of: " + min_anchor_mapq +  "for defining anchors");
		}
		
		
		PrintWriter outFq = null;
		SAMFileWriter outputSam = null;
		
		try {
			
			//TODO:
			SAMFileHeader samFileHeader = new BAMCollection(bams, samples).getMergedHeader(SAMFileHeader.SortOrder.unsorted);
			
			//TODO: Open up the .fq and .bam writer here, then run the potential MEIFinder
			outFq = new PrintWriter(new FileWriter(outfile + "_potential.fq"), true);
			outputSam = SAMWriting.makeSAMWriter(new File(outfile + "_potential.bam"), samFileHeader, new File(tmp), memory, SAMFileHeader.SortOrder.unsorted, true);
			
			//TODO: Loop over runPotentialMEIFinder depending on the number of BAM files given
			for (String file : bams){
				runPotentialMEIFinder(file, outFq, outputSam, tool, useSplit, minClipping, maxClipping);
			}
			
			

		} catch (IOException e) {
			logger.error("[PMRF] Could not create or find file / directory");
			logger.error(e.getMessage());
		} finally {
			if (outFq != null){
				outFq.close();
			}
			if (outputSam != null){
				outputSam.close();
			}
		}
		
		long end = System.currentTimeMillis();
		long millis = end - start;
		String time = String.format("%d min, %d sec", 
			    TimeUnit.MILLISECONDS.toMinutes(millis),
			    TimeUnit.MILLISECONDS.toSeconds(millis) - 
			    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
			);
		
		logger.info("[PMRF] PotentialMEIReadFinder ran in : " + time);
		
	}

	public static void runPotentialMEIFinder(String inFile, PrintWriter outFq, SAMFileWriter outSam, String mappingTool,
			boolean useSplit, int minClipping, int maxClipping) {
		
		File inBam = new File(inFile);
		PotentialMobilePairIterator potentialMEIReads = new PotentialMobilePairIterator(inBam, mappingTool, useSplit,
																minClipping, maxClipping, min_avg_qual, min_anchor_mapq);
		//SAMFileHeader samFileHeader = potentialMEIReads.getSAMReader().getFileHeader();
		
		//File tmpFile = new File(tmp);
		
		

		//If file already exists for 1st time then maybe program should issue an error?
		//PrintWriter outFq = new PrintWriter(new FileWriter(outPrefix + "_potential.fq"), true);
		//SAMFileWriter outputSam = SAMWriting.makeSAMWriter(new File(outPrefix + "_potential.bam"), samFileHeader, tmpFile, maxMemory, SAMFileHeader.SortOrder.unsorted, true);
		

		int c = 0;
		int d = 0;
	
		for (SAMRecordHolderPair<NrMappingsSAMRecordHolder> pair : YieldUtils.toIterable(potentialMEIReads)){

				pair.writeMobileReadsToFastQAndPotentialPairsToSAM(outFq, outSam, useSplit, true);
				c++;
				if(pair.hasSplitReadOfCertainSize()){
					d++;
				}

		}
	
		//outFq.close();
		//outputSam.close();
		logger.info("Found nr of potential mobile pairs supporting MEI events: " + c);
		logger.info("Of which " + d + " are of split-read nature.");


	}
	
	private static Options addCmdOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("BAM file containing mapped and unmapped paired-end reads" +
				" against human reference genome");
		
		options.addOption(OptionBuilder.create("in"));
		
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Output prefix for BAM files and fq files" +
				" containing supporting potential MEI paired-end reads");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("Tool Name");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Tool used for read mapping against reference");
		
		options.addOption(OptionBuilder.create("tool"));
		
		OptionBuilder.withDescription("use this option (-split) if you want to include split reads");
		
		options.addOption(OptionBuilder.create("split"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum size of clipped sequence for a split-read. Only used " +
				"when -split is supplied.");
		
		options.addOption(OptionBuilder.create("minclip"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum size of clipping for the end opposite of the end which has a " +
				"minimum clipping of at least the number specified by -minclip. Only used " +
				"when -split is supplied.");
		
		options.addOption(OptionBuilder.create("maxclip"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum average base quality for clipped read. Default = " + min_avg_qual);
		
		options.addOption(OptionBuilder.create("min_avg_qual"));
		
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
		OptionBuilder.withDescription("Minimum mapq needed for anchors. ONLY USED when -tool is unspecified. Default: " + min_anchor_mapq );
		
		options.addOption(OptionBuilder.create("mapq"));
		
		return options;
		
	}
}
