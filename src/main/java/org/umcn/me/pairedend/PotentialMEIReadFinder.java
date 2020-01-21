package org.umcn.me.pairedend;

import com.google.code.jyield.YieldUtils;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileWriter;
import org.apache.commons.cli.*;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.me.sam.PotentialMobilePairIterator;
import org.umcn.me.samexternal.*;
import org.umcn.me.util.BAMCollection;
import org.umcn.me.util.BAMSample;
import org.umcn.me.util.MobileDefinitions;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

public class PotentialMEIReadFinder {
	
	public static Logger logger = Logger.getLogger(PotentialMEIReadFinder.class.getName());
	private static int min_avg_qual = 20;
	private static int min_anchor_mapq = 20;
	private static boolean skip_um_pairs = false; //um --> unique - multiply mapped pairs
	private static boolean query_sort_input = false; //query sort input for lower RAM usage
	
	//TODO Command line calling now does not implement the multiple BAM input feature
	public static void main(String[] args) {
		
		BasicConfigurator.configure();
		
		Options options = addCmdOptions();

		if (args.length == 0) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("java -Xmx4g -jar PotentialMEIFinder.jar", options);
		} else {
			CommandLineParser parser = new GnuParser();
			try {

				long start = System.currentTimeMillis();
				CommandLine line = parser.parse(options, args);

				if (line.hasOption("properties")) {
					Properties props = new Properties();
					File propertiesFile = new File(line.getOptionValue("properties"));
					props.load(new FileInputStream(propertiesFile));
					ValuesParsedFromPropertiesOrCmdLine valuesParsedFromPropertiesFile = new ValuesParsedFromPropertiesOrCmdLine(props);
					valuesParsedFromPropertiesFile.log();
				} else {
					ValuesParsedFromPropertiesOrCmdLine valuesParsedFromCmdLine = new ValuesParsedFromPropertiesOrCmdLine(line);
					valuesParsedFromCmdLine.log();
				}
				
				//TODO: Open up the .fq and .bam writer here, then run the potential MEIFinder
				
				//Loop over runPotentialMEIFinder depending on the number of BAM files given
				//runPotentialMEIFinder(infile, outfile, tool, tmp, useSplit, minClipping, maxClipping, memory); 

				//TODO: Close the .fq and .bam writer here.
				
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

	private static Options addCmdOptions() {
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

	private static final class ValuesParsedFromPropertiesOrCmdLine {
		String infile;
		String outfile;
		String tool;
		File tmp;
		Boolean useSplit;
		int minClipping;
		int maxClipping;
		int memory;

		ValuesParsedFromPropertiesOrCmdLine (Properties props) throws IOException{
			infile = props.getProperty(MobileDefinitions.INFILE).trim();

			outfile = props.getProperty(MobileDefinitions.OUTFILE).trim();

			if (props.containsKey(MobileDefinitions.MAPPING_TOOL)) {
				tool = props.getProperty(MobileDefinitions.MAPPING_TOOL).trim();
			} else {
				//3-12-2014:Changed from mosaik to unspecified
				tool = SAMDefinitions.MAPPING_TOOL_UNSPECIFIED;
			}

			useSplit = Boolean.parseBoolean(props.getProperty(MobileDefinitions.USE_SPLIT).trim());

			if (props.containsKey(MobileDefinitions.TMP)) {
				tmp = new File(props.getProperty(MobileDefinitions.TMP).trim() + File.separator + "mob_" + Long.toString(System.nanoTime()));
			} else {
				tmp = new File(System.getProperty("java.io.tmpdir") + File.separator + "mob_" + Long.toString(System.nanoTime()));
			}
			if ( ! tmp.mkdir() ) {
				throw new IOException("Can not create tmp directory: " + tmp);
			}

			minClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_CLIPPING).trim());
			maxClipping = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_CLIPPING).trim());

			memory = Integer.parseInt(props.getProperty(MobileDefinitions.MEMORY).trim());

			min_avg_qual = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_QUAL).trim());

			if (props.containsKey(MobileDefinitions.GRIPS_DISCARD_UNIQUE_MULTIPLE)) {
				skip_um_pairs = Boolean.parseBoolean(props.getProperty(MobileDefinitions.GRIPS_DISCARD_UNIQUE_MULTIPLE).trim());
			}

			if (tool.equals(SAMDefinitions.MAPPING_TOOL_UNSPECIFIED)
					&&
					props.containsKey(MobileDefinitions.MIN_MAPQ_ANCHOR)) {
				min_anchor_mapq = Integer.parseInt(props.getProperty(MobileDefinitions.MIN_MAPQ_ANCHOR));
			}
		}

		ValuesParsedFromPropertiesOrCmdLine (CommandLine line) {
			infile = line.getOptionValue("in");
			outfile = line.getOptionValue("out");
			tool = line.getOptionValue("tool", SAMDefinitions.MAPPING_TOOL_UNSPECIFIED); //3-12-2014:Changed from mosaik to unspecified
			useSplit = line.hasOption("split");
			tmp = new File( line.getOptionValue("tmp", System.getProperty("java.io.tmpdir")) );
			minClipping = Integer.parseInt(line.getOptionValue("minclip", "35")); //30-1-2015: typo in option parsing min --> minclip
			maxClipping = Integer.parseInt(line.getOptionValue("maxclip", "7")); //30-1-2015: typo in option parsing max --> maxclip
			memory = Integer.parseInt(line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));
			min_avg_qual = Integer.parseInt(line.getOptionValue("min_avg_qual", Integer.toString(min_avg_qual)));
			min_anchor_mapq = Integer.parseInt(line.getOptionValue("mapq",Integer.toString(min_anchor_mapq)));

			//TODO: add in skip um pair as an cli option
		}

		void log() {
			logger.info("Using infile(s): " + infile);
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
		}
	}

	// =================================================================================================================

	//TODO remove duplicate code for this function
	public static void runFromProperties(Properties props) throws IOException {

		long start = System.currentTimeMillis();
		ValuesParsedFromPropertiesOrCmdLine valuesParsedFromPropertiesFile = new ValuesParsedFromPropertiesOrCmdLine(props);
		valuesParsedFromPropertiesFile.log();

		if (props.containsKey(MobileDefinitions.QUERY_SORT_INPUT)) {
			query_sort_input = Boolean.parseBoolean(props.getProperty(MobileDefinitions.QUERY_SORT_INPUT).trim());
		}

		//TODO: If multiple samples were found, turn on the multiple_sample_calling option in the properties from Mobster.java

		//TODO: Sample parsing
		String[] samples = props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0);
		String[] bams = props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0);

		//Check to see whether bam names are unique
		if (bams.length != new HashSet<String>(Arrays.asList(bams)).size()) {
			logger.error("Supplied bams are not unique");
			System.exit(1);
		}

		if (bams.length != samples.length) {
			logger.error("Number of sample names do not match number of supplied bams");
			System.exit(1);
		}

		logger.info("Number of samples detected: " + samples.length);

		PrintWriter outFq = null;
		SAMFileWriter outputSam = null;
		SAMFileHeader samFileHeader = null;
		try {
			
			outFq = new PrintWriter(new FileWriter(valuesParsedFromPropertiesFile.outfile + "_potential.fq"), true);
			
			//If more than one bam then try to make unique RG ids and associate sample name to each bam
			if (bams.length > 1) {
				logger.info("[PMRF] detected multiple samples, modifying RG's");
				BAMCollection collection = new BAMCollection(bams, samples, true);
				
				samFileHeader = collection.getMergedHeader(SAMFileHeader.SortOrder.unsorted);
				outputSam = SAMWriting.makeSAMWriter(new File(valuesParsedFromPropertiesFile.outfile + "_potential.bam"),
						samFileHeader, valuesParsedFromPropertiesFile.tmp, valuesParsedFromPropertiesFile.memory,
						SAMFileHeader.SortOrder.unsorted, true);
				
				//TODO: Loop over runPotentialMEIFinder depending on the number of BAM files given
				for (BAMSample bamSample : collection.getCloneOfBAMSampleList()){
					worker(bamSample.getBam(), valuesParsedFromPropertiesFile, outputSam, outFq,
							collection.getPrefixReadGroupIdFromBam(bamSample));
				}
			// Otherwise if there is just 1 bam, do not try to modify the read groups.
			} else if (bams.length == 1) {
				logger.info("[PMRF] detected one sample, not modifying RG's");
				SAMSilentReader reader = new SAMSilentReader(new File(bams[0]));
				samFileHeader = reader.getFileHeader();
				outputSam = SAMWriting.makeSAMWriter(new File(valuesParsedFromPropertiesFile.outfile + "_potential.bam"),
						samFileHeader, valuesParsedFromPropertiesFile.tmp, valuesParsedFromPropertiesFile.memory,
						SAMFileHeader.SortOrder.unsorted, true);
				reader.close();

				File bamSample = new File(bams[0]);
				worker(bamSample, valuesParsedFromPropertiesFile, outputSam, outFq, "");
			}
		} catch (IllegalArgumentException e) {
			logger.fatal(e.getMessage());
		} catch (IOException e) {
			logger.error("[PMRF] Could not create or find file / directory");
			logger.error(e.getMessage());
		} finally {
			if (outFq != null) {
				outFq.close();
			}
			if (outputSam != null) {
				outputSam.close();
			}
			if (valuesParsedFromPropertiesFile.tmp != null && ! valuesParsedFromPropertiesFile.tmp.delete() ) {
				logger.error("[PMRF] Could not delete temp: " + valuesParsedFromPropertiesFile.tmp);
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

	private static void worker(File bamSample, ValuesParsedFromPropertiesOrCmdLine values,
							   SAMFileWriter outputSam, PrintWriter outFq, String prefix) throws FileNotFoundException {
		//Query sort input if user wants this
		if (query_sort_input){
			File nameSortedBam = new File(bamSample.toString() + ".query_sorted");
			SAMWriting.writeSortedSAMorBAM(bamSample, nameSortedBam, values.tmp, values.memory, SortOrder.queryname);

			runPotentialMEIFinder(nameSortedBam.getAbsolutePath(), outFq, outputSam,
					values.tool, values.useSplit, values.minClipping, values.maxClipping, prefix);

			nameSortedBam.delete();
		}else{
			runPotentialMEIFinder(bamSample.getAbsolutePath(), outFq, outputSam,
					values.tool, values.useSplit, values.minClipping, values.maxClipping, prefix);
		}
	}

	private static void runPotentialMEIFinder(String inFile, PrintWriter outFq, SAMFileWriter outSam, String mappingTool,
			boolean useSplit, int minClipping, int maxClipping, String readGroupPrefix) {
		
		File inBam = new File(inFile);
		PotentialMobilePairIterator potentialMEIReads = new PotentialMobilePairIterator(inBam, mappingTool, useSplit,
																minClipping, maxClipping, min_avg_qual, min_anchor_mapq);
		logger.info("Skipping UM pairs?: " + skip_um_pairs);
		potentialMEIReads.setSkippingOfUMPairs(skip_um_pairs);

		//SAMFileHeader samFileHeader = potentialMEIReads.getSAMReader().getFileHeader();
		
		//File tmpFile = new File(tmp);

		//If file already exists for 1st time then maybe program should issue an error?
		//PrintWriter outFq = new PrintWriter(new FileWriter(outPrefix + "_potential.fq"), true);
		//SAMFileWriter outputSam = SAMWriting.makeSAMWriter(new File(outPrefix + "_potential.bam"), samFileHeader, tmpFile, maxMemory, SAMFileHeader.SortOrder.unsorted, true);

		int c = 0;
		int d = 0;
	
		for (SAMRecordHolderPair<NrMappingsSAMRecordHolder> pair : YieldUtils.toIterable(potentialMEIReads)){
				pair.writeMobileReadsToFastQAndPotentialPairsToSAM(outFq, outSam, useSplit, true, readGroupPrefix);
				c++;
				if (pair.hasSplitReadOfCertainSize()) {
					d++;
				}
            if ((c % 200000) == 0)
                logger.info("Processed " + c + " mobile read pairs");
            if ((d % 50000) == 0)
                logger.info("Processed " + d + " mobile read pairs that have split read or certain size");
		}
	
		//outFq.close();
		//outputSam.close();
		logger.info("Found nr of potential mobile pairs supporting MEI events: " + c);
		logger.info("Of which " + d + " are of split-read nature.");
	}
}
