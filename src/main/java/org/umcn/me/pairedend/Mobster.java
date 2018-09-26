package org.umcn.me.pairedend;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Properties;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.util.MobileDefinitions;

public class Mobster {

	
	public static Logger logger = Logger.getLogger("Mobster");
	public static String propertiesFile = null;
	public static String inFile = null;
	public static String outFile = null;
	public static String sampleName = null;
	public static final String VERSION;
	
	static{
		Properties prop = new Properties();
		try {
			prop.load(Mobster.class.getResourceAsStream("/properties/version.properties"));
		} catch (FileNotFoundException e) {
			logger.error("Could not load Mobster properties", e);
		} catch (IOException e) {
			logger.error("Could not read Mobster properties", e);
		}
		VERSION = prop.getProperty("version");
	}
	
	public static void main(String[] args) {
		
		checkParams(args);
		
		Properties props = new Properties();
		
		BasicConfigurator.configure();
		
		if (propertiesFile != null){
			try {
				
				props.load(new FileInputStream(propertiesFile));
				if (inFile != null){
					props.put(MobileDefinitions.INFILE, inFile);
				}
				if (outFile != null){
					props.put(MobileDefinitions.OUTFILE, outFile);
				}
				
				if (sampleName != null){
					props.put(MobileDefinitions.SAMPLE_NAME, sampleName);
				}
				
				if (props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0).length !=
						props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0).length){
					logger.fatal("Number of supplied samples does not equal the number of supplied bams. Exiting");
					System.exit(1);
				} else if (props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0).length > 1){
					logger.info("Detected multiple samples. Multiple sample calling will be turned ON even if this property was set to false in the properties file");
					props.put(MobileDefinitions.MULTIPLE_SAMPLE_CALLING, "true");
				}
				
				if (props.containsKey(MobileDefinitions.PICARD_COLLECT_INSERT_METRICS) &&
						"true".equals(props.getProperty(MobileDefinitions.USE_PICARD).trim())){
					
					//NOTE if multiple BAM Files are provided, only the insert size is investigated of the 1st BAM
					String picardCommand = "java -Xmx4g -jar " + props.getProperty(MobileDefinitions.PICARD_COLLECT_INSERT_METRICS) +
							" VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_hist.pdf" +
							" INPUT=" + props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP,0)[0] + " OUTPUT=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_insertstats" +
							" STOP_AFTER=50000000";

					execUnixCommand(picardCommand);
					
					//Values are on line 8
					//[0 = MEDIAN] [4 = MEAN] [17=99percentile]
					BufferedReader br = new BufferedReader(new FileReader(props.getProperty(MobileDefinitions.OUTFILE) + "_insertstats"));
					String line;
					boolean read = false;
					double median;
					int mean = 0;
					double percentile99 = 0;
					int clustermax = 0;
					int sd = 0;
					
					while ((line = br.readLine()) != null){
						if(read){
							String[] split = line.split("\t", -1);
							median = Double.parseDouble(split[0]);
							mean = (int) Double.parseDouble(split[4]);
							percentile99 = Double.parseDouble(split[17]);
							clustermax = (int) (median + percentile99);
							sd = (int) Double.parseDouble(split[5]);
							
							break;
						}else if(line.startsWith("MEDIAN_INSERT")){
							read = true;
						}
					}
					
					br.close();
					
					if (mean == 0 || clustermax == 0 || sd == 0){
						logger.error("Could not parse the PICARD CollectInsertSizeMetrics successfully");
						System.exit(-1);
					}
					
					if (props.containsKey(MobileDefinitions.USE_READ_LENGTH) && 
							"true".equalsIgnoreCase(props.getProperty(MobileDefinitions.USE_READ_LENGTH).trim())){
						clustermax = clustermax - Integer.parseInt(props.getProperty(MobileDefinitions.READ_LENGTH).trim());
					}
					
					props.put(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS, Integer.toString(clustermax));
					props.put(MobileDefinitions.MEAN_FRAGMENT_LENGTH, Integer.toString(mean));
					props.put(MobileDefinitions.SD_FRAGMENT_LENGTH, Integer.toString(sd));					
				}
				
				props.put(MobileDefinitions.INFILE_FROM_MOBIOME_MAPPING, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_mappedpotentials.bam");
				props.put(MobileDefinitions.INFILE_FROM_POTENTIAL_MEI_FINDER, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.bam");
				props.put(MobileDefinitions.ANCHOR_FILE, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_anchors.bam");
				props.put(MobileDefinitions.SPLIT_ANCHOR_FILE, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_splitanchors.bam");
				
				//do the mobiome mapping inhere
				String mobiomeMappingCmd = props.getProperty(MobileDefinitions.MOBIOME_MAPPING_CMD).trim();
				mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(FASTQ\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.fq");
				mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(DAT_FILE\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.dat");
				mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(OUT_FILE\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_mappedpotentials");

				PotentialMEIReadFinder.runFromProperties(props);
				
				int exitStatus = execUnixCommand(mobiomeMappingCmd);
				
				if (exitStatus != 0){
					logger.fatal("Mobiome mapping finished with non-normal exit value: " + exitStatus);
					logger.fatal("Therefore Mobster execution is terminated");
					System.exit(-1);
				}
				
				RefAndMEPairFinder.runFromPropertiesFile(props);
				
				AnchorClusterer.runFromPropertiesFile(props);
				
			} catch (IOException e) {
				logger.error(e.getMessage());
				System.exit(-1);
			}
		}
		
	}
	
	
	public static int execUnixCommand(String cmd){

		ProcessBuilder builder = new ProcessBuilder("/bin/sh", "-c", cmd);
		builder.redirectErrorStream(true);
		try {
			String s;
			Process process = builder.start();
            BufferedReader stdout = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
                while ((s = stdout.readLine()) != null) {
                    System.out.println(s);
                }
                int exitStatus = process.waitFor();
                
                System.out.println("Exit value: " + exitStatus);
                
                return exitStatus;
			
		} catch (IOException e) {
			logger.error("Error in executing command: " + cmd);
			logger.error(e.getMessage());
			return -1;
		} catch (InterruptedException e) {
			logger.error("Interrupted while executing command: " + cmd);
			logger.error(e.getMessage());
			return -2;
		}
	}
	
	protected static void checkParams(String[] args) {

		if (args.length > 0 && args[0] != null) {
			for (int i = 0; i < args.length; i = i + 2) {
				if (args[i].equalsIgnoreCase("-properties")){
					propertiesFile = args[i + 1];
				}else if (args[i].equalsIgnoreCase("-in")){
					inFile = args[i + 1];
				}else if (args[i].equalsIgnoreCase("-out")){
					outFile = args[i + 1];
				}else if (args[i].equalsIgnoreCase("-sn")){
					sampleName = args[i + 1];
				}
				else {
					logger.info("Invalid arguments. Please try again.");
					printUsage();
					throw new IllegalArgumentException("Invalid arguments. " + 
						"Acceptable arguments are one of [\"-properties\", \"-in\", \"-out\", \"-sn\"]");
				}
			}
		} else {
			logger.info("No argument provided. Please try again.");
			printUsage();
			throw new IllegalArgumentException("No argument provided.");
		}
	}
	
	public static void printUsage(){
		System.out.println("##########################");
		System.out.println("#MOBSTER                 #");
		System.out.println("##########################");
		System.out.println("Version: " + VERSION);
		System.out.println("Author: Djie Tjwan Thung");
		System.out.println("");
		System.out.println("Predict non-reference Mobile Element Insertion (MEI) events using one properties file.");
		System.out.println("\t-properties [properties]");
		System.out.println("\t-in [input .bam file]. This value will override corresponding value in properties file. Multiple BAM files may be specified if seperated by a comma");
		System.out.println("\t-out [output prefix]. This value will override corresponding value in properties file.");
		System.out.println(("\t-sn [sample name]. This value will override corresponding value in properties file. Multiple sample names may be specified if seperated by a comma"));
		System.out.println("Default mapping tool: " + SAMDefinitions.MAPPING_TOOL_UNSPECIFIED);
	}
}