package org.umcn.me.pairedend;

import java.io.IOException;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.umcn.gen.sam.SAMDefinitions;
import org.umcn.gen.sam.UnknownParamException;


public class PotentialMEIFinderMain {
	
	public static Logger logger = Logger.getLogger("PotentialMEIFinderMain");
	
	public static void main(String[] args){
		Options options;
		HelpFormatter formatter = new HelpFormatter();		
		BasicConfigurator.configure();
		PotentialMEIFinder potentialMobileFinder;
		
		options = addCmdOptions();
		String sep = "/";
		String infile;
		String outfile;
		String tool;
		
		if(args.length == 0){
			formatter.printHelp("java -Xmx4g -jar PotentialMEIFinder.jar", options);
		}else{
			CommandLineParser parser = new GnuParser();
			try {
				CommandLine line = parser.parse(options, args);
				if(line.hasOption("sep")){
					sep = line.getOptionValue("sep");
				}
				infile = line.getOptionValue("in");
				outfile = line.getOptionValue("out");
				tool = line.getOptionValue("tool", SAMDefinitions.MAPPING_TOOL_MOSAIK);
				
				logger.info("Running PotentialMEIFinder with follow settings:");
				logger.info("in: " + infile);
				logger.info("out: " + outfile);
				logger.info("sep: " + sep);
				logger.info("used mapping tool: " + tool);	
				
				potentialMobileFinder = new PotentialMEIFinder(sep, infile, outfile, tool);
				
				potentialMobileFinder.extractPotentialMobileReads();

				
			} catch (ParseException e) {
				logger.error("Error in parsing CLI arguments: " + e.getMessage());
			} catch (IOException e){
				logger.error("IO error: " + e.getMessage());
			} catch (UnknownParamException e){
				logger.error("Unknown mapping tool: " + e.getMessage());
			}
		}		
	}
	
	private static Options addCmdOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("BAM file containing mapped and unmapped paired-end reads" +
				" against human reference genome");
		
		options.addOption(OptionBuilder.create("in"));
		
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output prefix for BAM files and fq files" +
				" containing supporting potential MEI paired-end reads");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("read seperator");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Seperator used before read pair number 1 or 2. Default == '/'");
		
		options.addOption(OptionBuilder.create("sep"));
		
		OptionBuilder.withArgName("Tool Name");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Tool used for read mapping against reference");
		
		options.addOption(OptionBuilder.create("tool"));
		
		return options;
		
	}
}
