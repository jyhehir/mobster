package org.umcn.me.splitread;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.Logger;
import org.apache.log4j.BasicConfigurator;
import org.umcn.me.pairedend.RefAndMEPairFinder;
import org.umcn.me.sam.MobileSAMTag;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

/**
 * This class is used for extracting single-end reads from an original
 * bam file, where reads partially uniquely map to the reference
 * and partially uniquely map to the mobile reference.
 * 
 * Class is heavily dependent on methods of
 * org.umcn.me.pairedend.RefAndMEPairFinder.java
 * 
 * @author Djie
 *
 */
public class ExtractSplitAnchors {
	
	public static Logger logger = Logger.getLogger("ExtractSplitAnchors");
	
	public static void main(String[] args) {
		
		String bioBam;
		String mobBam;
		String outBam;
		String tool;
		
		Options options;
		
		BasicConfigurator.configure();
		HelpFormatter formatter = new HelpFormatter();
		options = createOptions();
		
		if (args.length == 0){
			formatter.printHelp("java -Xmx4g -jar ExtractSplitAnchors.jar", options);
		}else{
			CommandLineParser parser = new GnuParser();
			try {
				CommandLine line = parser.parse(options, args);
				bioBam = line.getOptionValue("biobam");
				mobBam = line.getOptionValue("mobbam");
				outBam = line.getOptionValue("out");
				tool = line.getOptionValue("tool");
				
				runExtractSplitAnchors(bioBam, mobBam, outBam, tool);
			} catch (ParseException e) {
				logger.fatal("Error when validating CLI arguments: " + e.getMessage());
			} catch (IOException e){
				logger.fatal("Error in reading / writing files: " + e.getMessage());
			}
		}
						
	}
	
	/**
	 * create CMD line options for this class
	 * @return
	 */
	private static Options createOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Lifescope / Bioscope bam file");
		
		options.addOption(OptionBuilder.create("biobam"));
		
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Bam file containing trimmed reads mapped against mobile ref. "+
		"Read ids should correspond between -biobam and -mobbam");
		
		options.addOption(OptionBuilder.create("mobbam"));
		
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output prefix for bam file containing split anchors");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Tool used for mobile mapping");
		
		options.addOption(OptionBuilder.create("tool"));
		
		return options;
	}
	
	/**
	 * Get all reads which mapped (with clipped end) to mobile reference
	 * and output them to a new BAM file
	 * @param bioBam original BAM file
	 * @param mobBam BAM file with reads mapping to mobile reference
	 * @param outBam output location file for all anchors
	 * @param tool which tool used for mobile mapping (mosaik or bwa)
	 * @throws IOException
	 */
	public static void runExtractSplitAnchors(String bioBam, String mobBam, String outBam, String tool) throws IOException{
		
		Map<String, MobileSAMTag> mobileInfoMap = new HashMap<String, MobileSAMTag>();
		
		mobileInfoMap = RefAndMEPairFinder.getReadsMappingToME(new File(mobBam), new Vector<String>(), tool);

		logger.info("Nr of split reads from getReadsMappingToME: " + mobileInfoMap.keySet().size());
		
		RefAndMEPairFinder.writeRefCoordinateFile(new File(bioBam), mobileInfoMap, outBam, false, "");
		
		
	}
}
