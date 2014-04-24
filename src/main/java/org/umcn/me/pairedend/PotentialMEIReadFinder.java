package org.umcn.me.pairedend;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
import org.umcn.gen.sam.NrMappingsSAMRecordHolder;
import org.umcn.gen.sam.SAMDefinitions;
import org.umcn.gen.sam.SAMRecordHolderPair;
import org.umcn.gen.sam.SAMWriting;
import org.umcn.me.sam.PotentialMobilePairIterator;
import org.umcn.me.util.MobileDefinitions;

import com.google.code.jyield.YieldUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;

public class PotentialMEIReadFinder {
	
	public static Logger logger = Logger.getLogger("PotentialMEIReadFinder");
	private static int min_avg_qual = 20;
	
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
						tool = SAMDefinitions.MAPPING_TOOL_MOSAIK;
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
					
					
				}else{
					infile = line.getOptionValue("in");
					outfile = line.getOptionValue("out");
					tool = line.getOptionValue("tool", SAMDefinitions.MAPPING_TOOL_MOSAIK);
					useSplit = line.hasOption("split");
					tmp = line.getOptionValue("tmp", System.getProperty("java.io.tmpdir"));
					minClipping = Integer.parseInt(line.getOptionValue("min", "35"));
					maxClipping = Integer.parseInt(line.getOptionValue("max", "7"));
					memory = Integer.parseInt(line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));
					min_avg_qual = Integer.parseInt(line.getOptionValue("min_avg_qual", Integer.toString(min_avg_qual)));
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
				
				
				runPotentialMEIFinder(infile, outfile, tool, tmp, useSplit, minClipping, maxClipping, memory); 
				
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
		
		
		String infile;
		String outfile;
		String tool;
		String tmp;
		Boolean useSplit;
		int minClipping;
		int maxClipping;
		int memory;
		
		long start = System.currentTimeMillis();
		
		infile = props.getProperty(MobileDefinitions.INFILE).trim();
		outfile = props.getProperty(MobileDefinitions.OUTFILE).trim();
		
		if (props.containsKey(MobileDefinitions.MAPPING_TOOL)){
			tool = props.getProperty(MobileDefinitions.MAPPING_TOOL).trim();
		}else{
			tool = SAMDefinitions.MAPPING_TOOL_MOSAIK;
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
		
		
		try {
			runPotentialMEIFinder(infile, outfile, tool, tmp, useSplit, minClipping, maxClipping, memory);
		} catch (FileNotFoundException e) {
			logger.error("[PMRF] Could not create or find file / directory");
			logger.error(e.getMessage());
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

	public static void runPotentialMEIFinder(String inFile, String outPrefix, String mappingTool, String tmp,
			boolean useSplit, int minClipping, int maxClipping, int maxMemory) throws FileNotFoundException {
		
		File inBam = new File(inFile);
		PotentialMobilePairIterator potentialMEIReads = new PotentialMobilePairIterator(inBam, mappingTool, useSplit,
																minClipping, maxClipping, min_avg_qual);
		SAMFileHeader samFileHeader = potentialMEIReads.getSAMReader().getFileHeader();
		
		File tmpFile = new File(tmp);
		
		
		try{
			PrintWriter outFq = new PrintWriter(new FileWriter(outPrefix + "_potential.fq"), true);
			SAMFileWriter outputSam = SAMWriting.makeSAMWriter(new File(outPrefix + "_potential.bam"), samFileHeader, tmpFile, maxMemory, SAMFileHeader.SortOrder.unsorted, true);
			

			int c = 0;
			int d = 0;
		
			for (SAMRecordHolderPair<NrMappingsSAMRecordHolder> pair : YieldUtils.toIterable(potentialMEIReads)){

					pair.writeMobileReadsToFastQAndPotentialPairsToSAM(outFq, outputSam, useSplit, true);
					c++;
					if(pair.hasSplitReadOfCertainSize()){
						d++;
					}

			}
		
			outFq.close();
			outputSam.close();
			logger.info("Found nr of potential mobile pairs supporting MEI events: " + c);
			logger.info("Of which " + d + " are of split-read nature.");

			}
		catch (IOException e){

		}
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
		
		return options;
		
	}
}
