package org.umcn.me.pairedend;

import com.sun.org.apache.xpath.internal.operations.Bool;
import net.sf.samtools.*;
import org.apache.commons.cli.*;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.me.output.Save;
import org.umcn.me.output.vcf.MobsterToVCF;
import org.umcn.me.sam.MateCluster;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.samexternal.SAMWriting;
import org.umcn.me.splitread.ClippedRead;
import org.umcn.me.splitread.InvalidHardClipCigarException;
import org.umcn.me.util.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * TODO There is now only one parameter (rpc = reads per cluster) which sets the minimum
 * nr of supporting reads for clusters belonging both to single and double clusters.
 * Preferably there should be a 2nd parameter which can send the threshold rpc for
 * single clusters differently than for double clusters.
 * 
 * move execCommand, createNSequence, writeDelimitedLine & writeDelimitedFile
 * to a more general util class
 * 
 * make a seperate container class for predictions
 */

/**
 * 
 * This class contains functionality to cluster reads that act as anchors
 * for mates that map to mobile elements (i.e. this is designed for paired-end reads).
 * 
 * IMPORTANT NOTE:
 * paired-end reads are assumed to be pointed inwards to each other, i.e.,
 * like with illumina technology :
 * 
 * ---------->
 * ATCGAATTATAGGCATCGGGGTATATCTAGGAGACCT
 *                             <--------
 *                            
 * Reads will be clustered if they are mapped on same strand and map to same mobile element super family.
 * 
 * Also prediction windows and insertion points are estimated.
 * This information is written to a prediction file.
 * 
 * 
 * #Version 0.1.6MULTI
 * Added multi sample support for bam files with read groups and sample names
 * 
 * @author Djie Thung
 *
 */
public class AnchorClusterer {
	
	public static Logger logger = Logger.getLogger(AnchorClusterer.class.getName());
	
	private static final int FILTER_REGION = 90;
	private static final String VERSION;
	
	private static int mean_frag_size = 470;
	private static int sd_frag_size = 35;
	private static int percentile_99_fragment = 500;
	private static int overlap_default = 50;
	private static int rpc_default = 2;
	private static int maxdist_default = 450;
	private static int clip_spacing = 15;
	private static int dup_max_size = 50;
	private static int del_max_size = 20;
	private static int search = 180;
	private static int min_total_hits = 5;
	private static int min_initial_cluster_size_split = 2;
	private static boolean multiple_sample_calling = false;
	private static boolean filter_by_read_counts_single_sample = false;
	private static boolean lenient_search = false;
	private static boolean grips_mode = false;
	private static String TMP;
	private static int memory;
	
	private static String repmask_file = "./hg19_alul1svaerv.txt";
	private static boolean usePrefixReference = true;
	
	private static Properties mobster_properties;
	
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
	
	public static Vector<MobilePrediction> runFromProperties(Properties props) throws IOException{
//		public static final String ANCHOR_FILE = "ANCHOR_BAM_FILE"; DONE
//		public static final String SPLIT_ANCHOR_FILE = "ANCHOR_SPLIT_BAM_FILE";
//		public static final String READS_PER_CLUSTER = "READS_PER_CLUSTER"; DONE
//		public static final String DISCORDANT_OVERLAP = "DISCORDANT_CLUSTERS_MAX_OVERLAP";
//		public static final String DISCORDANT_DISTANCE = "DISCORDANT_CLUSTER_MAX_DISTANCE";
//		public static final String MEAN_FRAGMENT_LENGTH = "MEAN_FRAGMENT_LENGTH";
//		public static final String SD_FRAGMENT_LENGTH = "SD_FRAGMENT_LENGTH";
//		public static final String LENGTH_99PROCENT_OF_FRAGMENTS = "LENGTH_99PROCENT_OF_FRAGMENTS";
//		public static final String MAX_SPACING_CLIPPED_READS = "MAX_SPACING_OF_CLIPPED_READS";
//		public static final String MAX_OVERLAP_CLIPPED_CLUSTERS = "MAX_OVERLAP_OF_CLIPPED_CLUSTERS";
//		public static final String MAX_DISTANCE_CLIPPED_CLUSTERS = "MAX_DISTANCE_OF_CLIPPED_CLUSTERS";
//		public static final String SEARCH_AREA = "NEIGHBORHOOD_WINDOW_BP";
//		public static final String MINIMUM_TOTAL_HITS = "MINIMUM_SUPPORTING_READS";
//		public static final String MINIMUM_INITIAL_SPLIT_CLUSTER_READS = "MINIMUM_INITIAL_SPLIT_CLUSTER_READS";
		
		//TMP
		//MEMORY
		
		File anchor;
		File anchorIndex;
		File clusterBam;
		File splitAnchor = null;
		File splitAnchorIndex = null;
		String sample = ""; // default empty
		String outPrefix;
		int rpc; //default 2
		int overlap; // default 50 
		int maxdist; // default 450

		//anchorIndex = new File(line.getOptionValue("in").replaceAll(".bam$", ".bai"));

		mobster_properties = props;
		
		anchor = new File(props.getProperty(MobileDefinitions.ANCHOR_FILE).trim());
		anchorIndex = new File(props.getProperty(MobileDefinitions.ANCHOR_FILE).trim().replaceAll(".bam$", ".bai"));
		clusterBam = new File(props.getProperty(MobileDefinitions.OUTFILE).trim() + "_clusters.bam");
		outPrefix = props.getProperty(MobileDefinitions.OUTFILE).trim();
		
		long start = System.currentTimeMillis();

		if (props.containsKey(MobileDefinitions.SAMPLE_NAME)){
			sample = props.getProperty(MobileDefinitions.SAMPLE_NAME).trim();
		}
		
		if (props.containsKey(MobileDefinitions.READS_PER_CLUSTER)){
			rpc = Integer.parseInt(props.getProperty(MobileDefinitions.READS_PER_CLUSTER).trim());
		}else{
			rpc = rpc_default;
		}
		
		if (props.containsKey(MobileDefinitions.DISCORDANT_OVERLAP)){
			overlap = Integer.parseInt(props.getProperty(MobileDefinitions.DISCORDANT_OVERLAP).trim());
		}else{
			overlap = overlap_default;
		}
		
		if (props.containsKey(MobileDefinitions.DISCORDANT_DISTANCE)){
			maxdist = Integer.parseInt(props.getProperty(MobileDefinitions.DISCORDANT_DISTANCE));
		}else{
			maxdist = maxdist_default;
		}
		
		if (props.containsKey(MobileDefinitions.MEAN_FRAGMENT_LENGTH)){
			mean_frag_size = Integer.parseInt(props.getProperty(MobileDefinitions.MEAN_FRAGMENT_LENGTH).trim());
		}
		
		if (props.containsKey(MobileDefinitions.SD_FRAGMENT_LENGTH)){
			sd_frag_size = Integer.parseInt(props.getProperty(MobileDefinitions.SD_FRAGMENT_LENGTH).trim());
		}
		
		if (props.containsKey(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS)){
			percentile_99_fragment = Integer.parseInt(props.getProperty(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.MAX_SPACING_CLIPPED_READS)){
			clip_spacing = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_SPACING_CLIPPED_READS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.MAX_DISTANCE_CLIPPED_CLUSTERS)){
			del_max_size = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_DISTANCE_CLIPPED_CLUSTERS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.MAX_OVERLAP_CLIPPED_CLUSTERS)){
			dup_max_size = Integer.parseInt(props.getProperty(MobileDefinitions.MAX_OVERLAP_CLIPPED_CLUSTERS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.SEARCH_AREA)){
			search = Integer.parseInt(props.getProperty(MobileDefinitions.SEARCH_AREA).trim());
		}
		
		if (props.containsKey(MobileDefinitions.MINIMUM_TOTAL_HITS)){
			min_total_hits = Integer.parseInt(props.getProperty(MobileDefinitions.MINIMUM_TOTAL_HITS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.MINIMUM_INITIAL_SPLIT_CLUSTER_READS)){
			min_initial_cluster_size_split = Integer.parseInt(props.getProperty(MobileDefinitions.MINIMUM_INITIAL_SPLIT_CLUSTER_READS).trim());
		}
		
		if (props.containsKey(MobileDefinitions.TMP)){
			TMP = props.getProperty(MobileDefinitions.TMP).trim() + File.separator + "mob_" + Long.toString(System.nanoTime());			
		}else{
			TMP = System.getProperty("java.io.tmpdir") + File.separator + "mob_" + Long.toString(System.nanoTime());
		}
		
		File tmp = new File(TMP);
		
		if ( ! tmp.mkdir() ){
			throw new IOException("Can not create tmp directory: " + tmp);
		}
		if (props.containsKey(MobileDefinitions.REPEAT_MASK_FILE)){
			repmask_file = props.getProperty(MobileDefinitions.RESOURCES_DIR) + props.getProperty(MobileDefinitions.REPEAT_MASK_FILE).trim();
		}
		if (props.containsKey(MobileDefinitions.PREFIX_REFERENCE)) {
			usePrefixReference = Boolean.getBoolean(props.getProperty(MobileDefinitions.PREFIX_REFERENCE));
		}
				
		if (props.containsKey(MobileDefinitions.SPLIT_ANCHOR_FILE)){
			splitAnchor = new File(props.getProperty(MobileDefinitions.SPLIT_ANCHOR_FILE).trim());
			splitAnchorIndex = new File(props.getProperty(MobileDefinitions.SPLIT_ANCHOR_FILE).trim().replaceAll(".bam$", ".bai"));
		}
		
		if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING) && "true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING).toLowerCase())){
			multiple_sample_calling = true;
		}
		if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT)  &&
				"true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT).toLowerCase())){
			filter_by_read_counts_single_sample = true;	
		}
		
		if (props.containsKey(MobileDefinitions.GRIPS_MODE) && "true".equals(props.getProperty(MobileDefinitions.GRIPS_MODE).toLowerCase())){
			grips_mode = true;
		}
		
		memory = Integer.parseInt(props.getProperty(MobileDefinitions.MEMORY).trim());

		logger.info("Using reads per cluster: " + rpc);
		logger.info("Using max overlap for double clusters: " + overlap);
		logger.info("Using max distance between 5 and 3 end clusters for double clusters of: " +
				maxdist);
		logger.info("using tmp: " + tmp);
		logger.info("using GRIPs mode: " + grips_mode);

		//Step 1: Clustering
		//Run simple clustering algorithm

		logger.info("Using simple clustering with search area of: " + search);
			
		runSimpleAnchorClusterer(anchor, anchorIndex, clusterBam, search, rpc, usePrefixReference);

		
		//Step 3: Make single & double cluster predictions out of the clusters
		Vector<MobilePrediction> matePredictions;
		
		matePredictions = clusterMateClusters(clusterBam, new File(clusterBam.toString().replaceAll(".bam$", ".bai")),
				overlap, maxdist);

		//Step 4: do split clustering when available
		if(splitAnchor != null){
			logger.info("Number of initial minimum hits for split clusters:" + min_initial_cluster_size_split);
			clusterSplitAnchorsToBAM(splitAnchor, splitAnchorIndex,
					new File(outPrefix + "_splitcluster.bam"), min_initial_cluster_size_split, usePrefixReference);

//			execCommand(samdir + "samtools sort " + line.getOptionValue("out") + "_splitcluster.bam" +
//			 " " + line.getOptionValue("out") + "_splitcluster_sorted");
//			execCommand(samdir + "samtools index " + line.getOptionValue("out") + "_splitcluster_sorted.bam");

			try {
				matePredictions = mergeSplitAndMateClusters(matePredictions, new File(outPrefix +
						"_splitcluster.bam"), new File(outPrefix + "_splitcluster.bai"));
			} catch (IllegalSAMPairException e) {
				logger.error(e.getMessage());
			}
		}

		matePredictions = mergePredictions(matePredictions, overlap, maxdist);

		long end = System.currentTimeMillis();
		long millis = end - start;
		String time = String.format("%d min, %d sec",
				TimeUnit.MILLISECONDS.toMinutes(millis),
				TimeUnit.MILLISECONDS.toSeconds(millis) -
				TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
			);

		logger.info("Anchorclusterer ran in : " + time);

		//Try to delete tmp
		if (tmp != null && ! tmp.delete() ){
			logger.error("Anchor Clusterer: Could not delete temp: " + tmp);
		}

		//Return predictions
		return matePredictions;
	}
	
	public static void main(String[] args) {
		Options options;
		Properties props = new Properties();
		
		BasicConfigurator.configure();
		HelpFormatter formatter = new HelpFormatter();
		options = createCmdOptions();
		
		if (args.length == 0){
			formatter.printHelp("java -Xmx4g -jar AnchorClusterer.jar", options);
		}else{
			CommandLineParser parser = new GnuParser();
			try {
				long start = System.currentTimeMillis();
				CommandLine line = parser.parse(options, args);

				//When a properties file has been supplied, load it
				if (line.hasOption("properties")){
					props.load(new FileInputStream(line.getOptionValue("properties")));
				}

				//Load the right command line arguments into the properties
				if(line.hasOption("insplit")) props.put(MobileDefinitions.SPLIT_ANCHOR_FILE, line.getOptionValue("insplit"));
				props.put(MobileDefinitions.SAMPLE_NAME, line.getOptionValue("sample", ""));
				props.put(MobileDefinitions.READS_PER_CLUSTER, line.getOptionValue("rpc", Integer.toString(rpc_default)));
				props.put(MobileDefinitions.DISCORDANT_OVERLAP, line.getOptionValue("overlap", Integer.toString(overlap_default)));
				props.put(MobileDefinitions.DISCORDANT_DISTANCE, line.getOptionValue("maxdist", Integer.toString(maxdist_default)));
				props.put(MobileDefinitions.MEAN_FRAGMENT_LENGTH, line.getOptionValue("mfl"));
				props.put(MobileDefinitions.SD_FRAGMENT_LENGTH, line.getOptionValue("sd"));
				props.put(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS, line.getOptionValue("maxclust", Integer.toString(percentile_99_fragment)));
				props.put(MobileDefinitions.SEARCH_AREA, line.getOptionValue("search", Integer.toString(search)));
				props.put(MobileDefinitions.MINIMUM_TOTAL_HITS, line.getOptionValue("mintotal", Integer.toString(min_total_hits)));
				props.put(MobileDefinitions.MINIMUM_INITIAL_SPLIT_CLUSTER_READS, line.getOptionValue("splithits", Integer.toString(min_initial_cluster_size_split)));
				props.put(MobileDefinitions.TMP, line.getOptionValue("tmp", System.getProperty("java.io.tmpdir")));
				props.put(MobileDefinitions.REPEAT_MASK_FILE, line.getOptionValue("repmask", repmask_file));
				props.put(MobileDefinitions.MULTIPLE_SAMPLE_CALLING, Boolean.toString(line.hasOption("multiplesample")));
				props.put(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT, line.hasOption("multisample_stringent"));
				props.put(MobileDefinitions.GRIPS_MODE, Boolean.toString(line.hasOption("grips")));
				props.put(MobileDefinitions.OUTFILE, line.getOptionValue("out") + "_predictions.txt");
				props.put(MobileDefinitions.VCF, Boolean.toString(line.hasOption("vcf")));
				props.put(MobileDefinitions.MEMORY, line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));

				props.put(MobileDefinitions.MAX_SPACING_CLIPPED_READS, Integer.toString(clip_spacing));
				props.put(MobileDefinitions.MAX_DISTANCE_CLIPPED_CLUSTERS, Integer.toString(del_max_size));
				props.put(MobileDefinitions.MAX_OVERLAP_CLIPPED_CLUSTERS, Integer.toString(dup_max_size));
				props.put(MobileDefinitions.PREFIX_REFERENCE, Boolean.toString(usePrefixReference));

				//Perform clustering
				Vector<MobilePrediction> predictions = runFromProperties(props);

				long end = System.currentTimeMillis();
				long millis = end - start;
				String time = String.format("%d min, %d sec", 
					    TimeUnit.MILLISECONDS.toMinutes(millis),
					    TimeUnit.MILLISECONDS.toSeconds(millis) -
					    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
					);
				
				logger.info("Anchorclusterer ran in : " + time);

				//Filter and save the result
				Save.runFromProperties(props, predictions);
				
			}  catch (IOException e){
				logger.error("Error in opening files: " + e.getMessage());
			} catch (org.apache.commons.cli.ParseException e) {
				// TODO Auto-generated catch block
				logger.error("Error in parsing CLI arguments: " + e.getMessage());
			}
			
		}

	}

	/**
	 * Executes a given argument literally on the command line
	 * @param cmd command to execute
	 */
	public static void execCommand(String cmd){
		Runtime rt = Runtime.getRuntime();
		try {
			Process p = rt.exec(cmd);
			int exitValue = p.waitFor();
			
			if (exitValue != 0){
				logger.error("Error in executing command (exit code = " + exitValue + " ): " + cmd);
				System.exit(-1);
			}
			p.getInputStream().close();
            p.getErrorStream().close();
            p.getOutputStream().close();
            p.destroy();
			
		} catch (IOException e) {
			logger.error("Error in executing command: " + cmd);
			logger.error(e.getMessage());
			System.exit(-2);
		} catch (InterruptedException e){
			logger.error("Error in executing command: " + cmd);
			logger.error(e.getMessage());
			System.exit(-3);
		}
		
	}
	
	public static Vector<MobilePrediction> mergePredictions(Vector<MobilePrediction> predictions, int overlapdef, int maxdist){
		Map<RegionWithLabel, MobilePrediction> regionMap = CollectionUtil.mobilePredictionsToRegions(predictions, overlap_default);
		List<RegionWithLabel> regions = new ArrayList<RegionWithLabel>(regionMap.keySet());
		Set<MobilePrediction> uniquePredictions = new HashSet<MobilePrediction>();
		
		Collections.sort(regions);
		int c = 0;

		for (RegionWithLabel r : regions){
			List<RegionWithLabel> overlap = r.overlapsWithSortedRegionList(regions);
			for (RegionWithLabel or : overlap){
				if (or.equals(r)){
					continue;
				}else{
					if ( mergeMobilePredictions(r, or, regionMap, overlapdef, maxdist)){
						c++;
					}
				}
			}
		}
		
		for (RegionWithLabel r : regions){
			uniquePredictions.add(regionMap.get(r));
		}
		
		logger.info("Number of merged Mobile Predictions: " + c);
		logger.info("Number of predictions before merging: " + predictions.size());
		logger.info("Number of predictions after merging: " + uniquePredictions.size());
		
		
		Map<RegionWithLabel, MobilePrediction> regionMap2 = CollectionUtil.mobilePredictionsToRegions(new Vector<MobilePrediction>(uniquePredictions), overlap_default);
		List<RegionWithLabel> regions2 = new ArrayList<RegionWithLabel>(regionMap2.keySet());
		Collections.sort(regions2);
		Vector<MobilePrediction> mergedPredictions =  new Vector<MobilePrediction>();
		for (RegionWithLabel r : regions2){
			mergedPredictions.add(regionMap2.get(r));
		}
		
		return mergedPredictions;
		
	}
	
	private static boolean mergeMobilePredictions(RegionWithLabel r1, RegionWithLabel r2,
			Map<RegionWithLabel, MobilePrediction> labeledToMobile, int overlapdef, int maxdist){
		
		MobilePrediction m1 = labeledToMobile.get(r1);
		MobilePrediction m2 = labeledToMobile.get(r2);
		
		Set<String> m1SplitIds = m1.getSplitClusterIds();		
		Set<String> m2SplitIds = m2.getSplitClusterIds();
		
		//We can not merge MobilePredictions containing different split reads
		if (m1SplitIds.size() >= 1 && m2SplitIds.size() >= 1 && ! m1SplitIds.equals(m2SplitIds)){
			return false;
		}
		
		boolean success = false;
		try {
			success = m1.addDiscordantReadsOfMobilePrediction(m2, overlapdef, maxdist);
		} catch (IllegalSAMPairException e) {
			logger.warn("Illegal SAM information when merging: " + e.getMessage());
		}
		if (success){
			labeledToMobile.put(r2, m1);
		}
		
		return success;

	}
	
	public static Options createCmdOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("String");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Properties file");
		
		options.addOption(OptionBuilder.create("properties"));
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("BAM file (COORDINATE SORTED!) containing mapped anchors");
		
		options.addOption(OptionBuilder.create("in"));
		
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Prefix for output files");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Mean fragment length");
		
		options.addOption(OptionBuilder.create("mfl"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Standard deviation of fragment length distribution");
		
		options.addOption(OptionBuilder.create("sd"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Reads per cluster: Minimum number of reads a cluster must have" +
				" (default is " + rpc_default + ")");
		
		options.addOption(OptionBuilder.create("rpc"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Number of total required reads for a mobile prediction (default: "+
				min_total_hits + ")");
		
		options.addOption(OptionBuilder.create("mintotal"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Search region for simple clustering algorithm. Anchors within the search" +
				" region will be added to the same cluster. (default is " + search + " )");
		
		options.addOption(OptionBuilder.create("search"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum distance between 5 end and 3 end cluster for combining them into" +
				" a double cluster (default is " + maxdist_default + "(bp))");
		
		options.addOption(OptionBuilder.create("maxdist"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Maximum overlap between 5 end and 3 end cluster for combining them into" +
				" a double cluster (default is " + overlap_default + " (bp)");
		
		options.addOption(OptionBuilder.create("overlap"));
		
//		OptionBuilder.withArgName("dir");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Directory to samtools. Default is '' (use this if samtools is in your PATH.");
//		
//		options.addOption(OptionBuilder.create("samdir"));
		
		OptionBuilder.withArgName(".out");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("RepeatMasker .out file for MEI filtering, default location: " +
				repmask_file);
		
		options.addOption(OptionBuilder.create("repmask"));
		
		OptionBuilder.withArgName("Boolean");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("If chromosomes should have the prefix 'chr' added to their dictionary " +
				usePrefixReference);
		
		options.addOption(OptionBuilder.create("usePrefixReference"));
		
		
		OptionBuilder.withArgName(".bam");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Coordinate sorted split anchor file");
		
		options.addOption(OptionBuilder.create("insplit"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Min number of initial split cluster hits");
		
		options.addOption(OptionBuilder.create("splithits"));
		
//		OptionBuilder.withDescription("Switch for doing lenient search. Meaning: " +
//				"Single cluster split reads are allowed," +
//				" with a minimum of " + min_total_hits +
//				". Single cluster split reads are not allowed to map against poly A or poly T." +
//				"Default use of lenient search: " + lenient_search);
//		
//		options.addOption(OptionBuilder.create("lenient"));
		
		OptionBuilder.withArgName("Integer");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Max expected cluster size. Preferably based on 99th width percentile outputted " +
				"by Picard's CollectInsertSizeMetrics");
		
		options.addOption(OptionBuilder.create("maxclust"));
		
		OptionBuilder.withArgName("String");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Sample name");
		
		options.addOption(OptionBuilder.create("sample"));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Tmp directory to use, when writing large files");
		
		options.addOption(OptionBuilder.create("tmp"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Alter the number of records in memory (increase to decrease number of file handles). Default: " + SAMWriting.MAX_RECORDS_IN_RAM);
		
		options.addOption(OptionBuilder.create("max_memory"));
		
		OptionBuilder.withDescription("Turn this option on if you want to do multiple sample calling. Note that the bam file should have @RG tags and each @RG tag should have a SM (sample field)");
		
		options.addOption(OptionBuilder.create("multiplesample"));
		
		OptionBuilder.withDescription("Only checked if -multiplesample is used. Filter predictions by read counts from a " +
		"single sample. I.e. if a prediction is supported by 3 reads from sample1 and 4 reads by sample2 and -mintotal equals 5, " +
				"then the prediction is filtered because both samples do not have 5 supporting reads. If the option is not turned on, " +
				"the prediction is kept because the 3 reads from sample1 and 4 reads from sample2 are added.");
		
		options.addOption(OptionBuilder.create("multisample_stringent"));
		
		
		Option grip = new Option("grips", "Change settings specifically for GRIP clustering");
		options.addOption(grip);

		Option vcf = new Option("vcf", "Output the results of Mobster in VCF format in addition to the default format.");
		options.addOption(vcf);

		return options;
	}
	
	/**
	 * Clustering using simple algorithm ;-)
	 * 
	 * This clustering method clusters anchors into clusters. Reads
	 * within the user-specified searchRegion, which have the same strand mapping
	 * and support the same mobile element family are clustered together. When a cluster
	 * has at least minReads number or reads in it and no more reads are within the searchRegion
	 * then a cluster is definite and will be written to a new bam file, which has information
	 * about the number of reads supporting the cluster, the mobile element family the cluster
	 * supports and the cluster length.
	 * 
	 * @param anchor BAM file containing anchor reads for mates mapping to mobile ref.
	 * Should be coordinate sorted
	 * @param index index file associated with anchor BAM file.
	 * @param clusterBam File to write found clusters to in BAM format.
	 * @param searchRegion all reads within this search region will be added to cluster,
	 * the reads do have to have same strand mapping and support same mobile element family.
	 * @param minReads number of reads a cluster has to have, otherwise potential cluster is
	 * thrown away.
	 */
	public static void runSimpleAnchorClusterer(File anchor, File index, File clusterBam, int searchRegion, int minReads, boolean prepend){

		SAMFileReader input = new SAMSilentReader(anchor, index);
		SAMFileHeader header = input.getFileHeader();
		SampleBam sampleCalling = SampleBam.SINGLESAMPLE;
		
		header = editSAMFileHeader(header, prepend);
		SAMFileWriter outputSam = SAMWriting.makeSAMWriter(clusterBam, header, new File(TMP), memory, SAMFileHeader.SortOrder.coordinate);
		
		Vector<MateCluster<SAMRecord>> mobileClusters = new Vector<MateCluster<SAMRecord>>();
		
		int c = 0;
		int skippedClustersBecauseOfNotSameRefMapping = 0;
		double minPercentSameMateRefMapping = 0.0;
		int maxDiffMateMapping = Integer.MAX_VALUE;
		boolean usePrefixReference = true;
		
		if (mobster_properties.containsKey(MobileDefinitions.GRIPS_MIN_PERCENT_SAME_REF_MATE_MAPPING)){
			minPercentSameMateRefMapping = Double.parseDouble(mobster_properties.getProperty(MobileDefinitions.GRIPS_MIN_PERCENT_SAME_REF_MATE_MAPPING));
		}
		
		if (mobster_properties.containsKey(MobileDefinitions.GRIPS_MAX_DIFF_MATE_MAPPINGS)){
			maxDiffMateMapping = Integer.parseInt(mobster_properties.getProperty(MobileDefinitions.GRIPS_MAX_DIFF_MATE_MAPPINGS));
		}
		
		if (mobster_properties.containsKey(MobileDefinitions.PREFIX_REFERENCE))
			usePrefixReference = Boolean.getBoolean(mobster_properties.getProperty(MobileDefinitions.PREFIX_REFERENCE));
		
		logger.info("Using a minPercentSameMateRefMapping threshold of: " + minPercentSameMateRefMapping);
		logger.info("Using a ZZdifferent mate chr threshold of: " + maxDiffMateMapping);
		
		if (multiple_sample_calling){
			sampleCalling = SampleBam.MULTISAMPLE;
		}
		
		mobileClusters.add(new MateCluster<SAMRecord>(header, false, true, sampleCalling));
		Iterator<MateCluster<SAMRecord>> iter = mobileClusters.iterator();
		for(SAMRecord record : input){
			boolean added = false;
			while (iter.hasNext()){
				MateCluster<SAMRecord> currentCluster = iter.next();
				if(currentCluster.isWithinSearchArea(record, searchRegion)){
					if(currentCluster.add(record)){
						added = true;
						break;
					}
				}else if(currentCluster.size() >= minReads){
					if (currentCluster.getHighestPercentageOfMateAlignmentsToSameChrosome(true) >= minPercentSameMateRefMapping
							&& currentCluster.getNumberOfDifferentChromosomeMappingsOfMates(true) <= maxDiffMateMapping){
						currentCluster.writeClusterToSAMWriter(outputSam, Integer.toString(c), usePrefixReference);
						c++;
					}else{
						skippedClustersBecauseOfNotSameRefMapping++;
					}
					iter.remove();
				}else{
					iter.remove();
				}
			}
			
			if (!added){
				MateCluster<SAMRecord> newCluster = new MateCluster<SAMRecord>(header, false, true, sampleCalling);
				newCluster.add(record);
				mobileClusters.add(newCluster);
			}
			iter = mobileClusters.iterator();
		}
		

		//Write anchors to clusters.bam which have not been written yet.
		for (MateCluster<SAMRecord> cluster : mobileClusters){
			if (cluster.size() >= minReads){
				cluster.writeClusterToSAMWriter(outputSam, Integer.toString(c), usePrefixReference);
			}
		}
		logger.info("Number of skipped clusters because mates do not align to same chromosome:" + skippedClustersBecauseOfNotSameRefMapping);
		logger.info("Found nr of clusters: " + c);

		input.close();
		outputSam.close();
	}
	
	public static void clusterSplitAnchorsToBAM(File anchor, File index, File outBam, int minClustersize, boolean prepend){
		
		SampleBam sampleCalling = SampleBam.SINGLESAMPLE;
		
		if (multiple_sample_calling){
			sampleCalling = SampleBam.MULTISAMPLE;
		}
		
		Set<String> alreadyClustered = new HashSet<String>();
		ClippedReadSet<ClippedRead> cluster = new ClippedReadSet<ClippedRead>(sampleCalling);
		
		SAMFileReader input = new SAMSilentReader(anchor, index);
		SAMFileReader input2 = new SAMSilentReader(anchor, index);
		
		SAMFileHeader header = input.getFileHeader();
		SAMFileHeader header2 = input2.getFileHeader();
		
		int c = 0;
		int d = 0;
		
		//convert reference names in SAM header from 1 to -> chr1 etc.
		header = editSAMFileHeader(header, prepend);
		header2 = editSAMFileHeader(header2, prepend);
		SAMFileWriter outputSam = SAMWriting.makeSAMWriter(outBam, header, new File(TMP), memory, SAMFileHeader.SortOrder.coordinate);
				
		for (SAMRecord record : input){
			
			try {
				ClippedRead read = new ClippedRead(record, true, false);
				c++;
				cluster = clusterClippedReads(alreadyClustered, input2, read);
				
				if ((cluster.size() - cluster.getNrOfLeftClippedPolyAMapping() - cluster.getNrOfLeftClippedPolyTMapping()
						- cluster.getNrOfRightClippedPolyAMapping() - cluster.getNrOfRightClippedPolyTMapping() ) >= minClustersize){
					d += cluster.size();
					cluster.writeClippedReadSetToBAM(outputSam);
				}
			} catch (InvalidHardClipCigarException e) {
				logger.warn("Invalid clipped read encountered when clustering split anchors: " + e.getMessage());
			}
		}
		
		logger.info("Nr of split reads found: " + c);
		logger.info("Nr of split reads clustered: " + d);
		outputSam.close();
		input.close();
		input2.close();
	}

	private static ClippedReadSet<ClippedRead> clusterClippedReads(Set<String> alreadyClustered,
			SAMFileReader input2, ClippedRead read)
			throws InvalidHardClipCigarException {
		
		SampleBam sampleCalling = SampleBam.SINGLESAMPLE;
		
		if (multiple_sample_calling){
			sampleCalling = SampleBam.MULTISAMPLE;
		}
		
		ClippedReadSet<ClippedRead> splitCluster = new ClippedReadSet<ClippedRead>(sampleCalling);
		SAMRecordIterator iter;
		if (!alreadyClustered.contains(read.getSAMRecord().getReadName())) {
			int clippedEnd = read.getClippedPosition();
			String ref = read.getSAMRecord().getReferenceName();
			splitCluster.add(read);
			alreadyClustered.add(read.getSAMRecord().getReadName());
			if (!read.isStartClipped()) {
				iter = input2.query(ref, clippedEnd
						- dup_max_size, clippedEnd + del_max_size, false);
			}else{
				iter = input2.query(ref, clippedEnd - del_max_size,
						clippedEnd + dup_max_size, false);
			}
				while (iter.hasNext()) {
					ClippedRead readToBeAdded = new ClippedRead(iter.next(),
							true, false);
					String readNameToBeAdded = readToBeAdded.getSAMRecord()
							.getReadName();

					if (alreadyClustered.contains(readNameToBeAdded)) {
						continue;
					}

					if (readToBeAdded.isStartClipped()) {
						if (splitCluster.isOnSameReference(readToBeAdded)
								&& (splitCluster
										.isSamePrimaryMobileMapping(readToBeAdded, false) ||
										splitCluster.hasSameHomoPolymerMapping(readToBeAdded))
								&& splitCluster.isLeftClippedReadWithinRegion(
										readToBeAdded, clip_spacing)) {
							alreadyClustered.add(readNameToBeAdded);
							splitCluster.add(readToBeAdded);
						}
					} else if (splitCluster.isOnSameReference(readToBeAdded)
							&& (splitCluster
									.isSamePrimaryMobileMapping(readToBeAdded, false) ||
									splitCluster.hasSameHomoPolymerMapping(readToBeAdded))
							&& splitCluster.isRightClippedReadWithinRegion(
									readToBeAdded, clip_spacing)) {
						splitCluster.add(readToBeAdded);
						alreadyClustered.add(readNameToBeAdded);
					}
				}
				iter.close();
		}
		return splitCluster;
	}
	
	
	
	private static SAMFileHeader editSAMFileHeader(SAMFileHeader header, boolean prepend) {
		List<SAMSequenceRecord> sequenceList = new ArrayList<SAMSequenceRecord>();
		
		for ( SAMSequenceRecord sequenceRec : header.getSequenceDictionary().getSequences()){
			StringBuilder sequenceBuilder = new StringBuilder();
			SAMSequenceRecord newSeq;
			if(!sequenceRec.getSequenceName().startsWith("chr") && prepend){
				sequenceBuilder.append("chr");
			}
			sequenceBuilder.append(sequenceRec.getSequenceName());
			newSeq = new SAMSequenceRecord(sequenceBuilder.toString(), sequenceRec.getSequenceLength());
			sequenceList.add(newSeq);
		}
		header.setSequenceDictionary(new SAMSequenceDictionary(sequenceList));
		
		return header;
	}

	
	/**
	 * Function which calls all other function after using a certain clustering method like
	 * runSimpleAnchorClusterer(). The BAM file containing clusters should be coordinate
	 * sorted and indexed.
	 * This function determines whether clusters can be formed together to form double clusters
	 * or whether they stay single clusters.
	 * Double ended clusters are made when there is one cluster on the forward strand and
	 *  one cluster on the reverse strand and they both support the same mobile element family.
	 *  Forward strand cluster should have lower genomic coordinate than reverse strand cluster.
	 *  Reverse strand cluster should be maximum of maxdist xx basepairs away
	 *   (this maxdist can be specified by CLI argument) and may overlap for xx bp (overlap arg)
	 *   at the end of the forward strand cluster. 
	 * 
	 * After making single and double clusters predictions are made and written to a file.
	 * 
	 * @param clusterIn
	 * @param clusterIndex
	 * @param overlap integer determining how much the inner border of a 3' end cluster may
	 * overlap with the inner border of a 5' end cluster.
	 * @param maxdist maximum distance between forward and reverse clusters to still cluster them as a
	 * double cluster.
tyuasx	 * @throws IOException
	 */
	public static Vector<MobilePrediction> clusterMateClusters(File clusterIn, File clusterIndex, int overlap, int maxdist) {
		
		SAMFileReader input = new SAMSilentReader(clusterIn, clusterIndex);
		SAMFileReader input2 = new SAMSilentReader(clusterIn, clusterIndex);	
		
		HashSet<String> doubleClusterNames = new HashSet<String>();
		Vector<MobilePrediction> mobilePredictions = new Vector<MobilePrediction>(); 
		
		int d = 0;
		int c = 0;
		int e = 0;
		for(SAMRecord cluster : input){
			try {
				MobilePrediction pred = makeMobilePredictionFromMateClusters(cluster, input2, maxdist, overlap, doubleClusterNames);
				if (pred != null && pred.getRightTotalHits() + pred.getLeftTotalHits() >= 2){
					mobilePredictions.add(pred);
					e++;
					if (pred.hasLeftMateCluster() && pred.hasRightMateCluster()){
						c++;
					}else{
						d++;
					}
				}
			} catch (IllegalSAMPairException ex) {
				// should not occur
			}			
		}

		logger.info("Nr of found double clusters: " + c);
		logger.info("Nr of found single clusters: " + d);
		logger.info("Nr of found clusters total: " + e);

		//Make prediction estimates and filter them when they overlap with annotated mobile elements
//		predictions = makePredictionVector(clusters);
//		filtered = filterKnownMEs(getKnownMEs(), predictions);
//		
//		logger.info("Nr of filtered predictions: " + (predictions.size() - filtered.size()));
//		
//		//Write all filtered predictions to file
//		String[] headerArray = {"Chr", "Mobile Element", "Insert Point", "border5", "border3",
//		"cluster5 length", "cluster3 length", "cluster5 hits", "cluster3 hits"};
//
//		Vector<String> head = new Vector<String>(Arrays.asList(headerArray));
//
//		writeDelimitedFile(filtered, head, predictionFile, "\t");	

		input.close();
		input2.close();
		
		return mobilePredictions;
	}
	
	public static MobilePrediction makeMobilePredictionFromMateClusters(SAMRecord cluster, SAMFileReader input,
			int maxdist, int overlap, HashSet<String> alreadyClustered) throws IllegalSAMPairException{
		
		Boolean reverse;
		SAMRecordIterator iter;
		SAMRecord clusterMate;
		String mobileFamily;
		String mobileFamilyMate;
		String ref;
		MobilePrediction prediction = null;
		int end;
		
		reverse = cluster.getReadNegativeStrandFlag();
		mobileFamily = cluster.getAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT).toString();
		ref = cluster.getReferenceName();
		end = cluster.getAlignmentEnd();

		if(!reverse){
			iter = input.query(ref, end - overlap, end + maxdist, false);
			while(iter.hasNext()){
				clusterMate = iter.next();
				mobileFamilyMate = clusterMate.getAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT).toString();
				if(clusterMate.getReadNegativeStrandFlag() && mobileFamily.equals(mobileFamilyMate)
						&& clusterMate.getAlignmentStart() >= end - overlap){

					prediction = new MobilePrediction(mean_frag_size, sd_frag_size, percentile_99_fragment, cluster, clusterMate);
					alreadyClustered.add(clusterMate.getReadName());
					alreadyClustered.add(cluster.getReadName());
	
					break;
				}
			}
			iter.close();
		}
		
		//If no double cluster was found add the 5'end or 3'end single cluster
		if(!alreadyClustered.contains(cluster.getReadName())){
			prediction = new MobilePrediction(mean_frag_size, sd_frag_size, percentile_99_fragment, cluster);
		}
		
		return prediction;
	}
	
	public static Vector<MobilePrediction> mergeSplitAndMateClusters(Vector<MobilePrediction> matePredictions,
			File splitClusterBam, File splitclusterIndex) throws IllegalSAMPairException{
		
		int mergedCount = 0;
		int multipleSplitClusters = 0;
		
		Set<String> alreadyClusteredReadNames = new HashSet<String>();
		SAMFileReader splitReader = new SAMSilentReader(splitClusterBam, splitclusterIndex);
				
		for (MobilePrediction matePrediction : matePredictions){
			String mobileHit = matePrediction.getMateMobileMapping();
			int left = matePrediction.getLeftPredictionBorder(50);
			int right = matePrediction.getRightPredictionBorder(50);
			
			SAMRecordIterator iter = splitReader.query(matePrediction.getOriginalReference(), left, right, false);
			MobilePrediction bestPrediction = null;
			int splitClustersInArea = 0;
			while (iter.hasNext()){
				splitClustersInArea++;
				MobilePrediction currentPrediction = new MobilePrediction(0, 0, 0, iter.next());
				if(currentPrediction.getMobileMappings().contains(mobileHit)){
					bestPrediction = highestConfidenceSplitPrediction(bestPrediction, currentPrediction);
				}
			}
			iter.close();
			if (splitClustersInArea >= 2){
				logger.info("This mate cluster region has more than 2 split clusters: " +
						matePrediction.getOriginalReference() +	":"+left+"-"+right);
				multipleSplitClusters++;
			}
			if (bestPrediction != null){
				mergedCount++;
				matePrediction.parseSAMRecordCluster(bestPrediction.getSplitSAMRecord());
				alreadyClusteredReadNames.add(bestPrediction.getSplitSAMRecord().getReadName());
			}
		}		
		splitReader.close();
		logger.info("Total mate predictions which have more than 1 split clusters in the vicinity: " + multipleSplitClusters);
		logger.info("Nr of mate clusters with merged split clusters: " + mergedCount);
		
		SAMFileReader splitReader2 = new SAMSilentReader(splitClusterBam, splitclusterIndex);
		
		int soloSplitClusters = 0;
		for (SAMRecord rec : splitReader2){
			if (!alreadyClusteredReadNames.contains(rec.getReadName())){
				MobilePrediction pred = new MobilePrediction(0, 0, 0, rec);
				if(pred.hasLeftAlignedSplitCluster() && pred.hasRightAlignedSplitCluster() && pred.hasOneOrZeroHomoPolymerMappings()){
					matePredictions.add(pred);
					soloSplitClusters++;
				}else if(lenient_search && pred.getLeftTotalHits() + pred.getRightTotalHits() >= min_total_hits
						&& pred.getLeftPolyAhits() == 0 && pred.getLeftPolyThits() == 0  && pred.getRightPolyAhits() == 0 
						&& pred.getRightPolyThits() == 0 ){
					matePredictions.add(pred);
					soloSplitClusters++;
				}
			}
		}
		
		logger.info("Total nr of split clusters without mate support: " + soloSplitClusters);
		splitReader2.close();
		//TODO filter predictions which do not meet requirements
		
		return matePredictions;
	}

	private static MobilePrediction highestConfidenceSplitPrediction(
			MobilePrediction prediction1, MobilePrediction prediction2) {
		//TODO if total hits cluster 1 == total hits cluster 2; assign highest confidence
		//split prediction to the one who is closer to the insertion estimate of mate cluster
		if(prediction1 == null){
			return prediction2;
		}else if(prediction2 == null){
			return prediction1;
		}
		
		boolean doubleCluster1 = prediction1.hasLeftAlignedSplitCluster() && prediction2.hasRightAlignedSplitCluster();
		boolean doubleCluster2 = prediction2.hasLeftAlignedSplitCluster() && prediction2.hasRightAlignedSplitCluster();
		int hitsCluster1 = prediction1.getLeftTotalHits() + prediction1.getRightTotalHits();
		int hitsCluster2 = prediction2.getLeftTotalHits() + prediction2.getRightTotalHits();
		
		if(doubleCluster1 && doubleCluster2){
			if(hitsCluster1 >= hitsCluster2){
				return prediction1;
			}else{
				return prediction2;
			}
		}else if(doubleCluster1){
			return prediction1;
		}else if(doubleCluster2){
			return prediction2;
		}else if(hitsCluster1 >= hitsCluster2){
			return prediction1;
		}else{
			return prediction2;
		}
		
		
	}
	
	/**
	 * This function filters predictions when the prediction window overlaps with an annotated
	 * mobile element stored in the parameter knownMEs.
	 * @param knownMEs - Chromosomes as key
     *                 - Dictionary as value: with mobile element super family as key
	 *                        and position start (0-based inclusive) and position end (0-based exclusive)
	 *                        as sting value. Start and end are seperated by -
	 * @param predictions Vector containing vector of predictions. Each prediction is a string with on
	 * index 0: reference chromosome for MEI
	 * index 1: mobile family for MEI
	 * index 2: estimated insert point
	 * index 3: border5 prediction window
	 * index 4: border3 prediction window
	 * index 5: cluster length of forward cluster (NA if not present)
	 * index 6: cluster length of reverse cluster (NA if not present)
	 * index 7: number of supporting reads for forward cluster (NA if not present)
	 * index 8: number of supporting reads for reverse cluster (NA if not present)
	 * @return prediction vector of the same format as arg predictions. Only now with predictions
	 * that do not overlap with elements kept in argument knownMEs.
	 */
	
	public static Vector<MobilePrediction> filterKnownMEs(Map<String, HashMap<String, HashSet<String>>> knownMEs,
			Vector<MobilePrediction> predictions){
		String ref;
		//String me;
		int border5;
		int border3;
		int alu = 0;
		int sva = 0;
		int l1 = 0;
		int herv = 0;
		int startKnown;
		int endKnown;
		boolean inKnownMe = false;
		
		Vector<MobilePrediction> filteredPredictions = new Vector<MobilePrediction>();
		
		for(MobilePrediction prediction : predictions){
			ref = prediction.getOriginalReference();
			//me = prediction.getMateMobileMapping();
			border5 = prediction.getLeftPredictionBorder(FILTER_REGION);
			border3 = prediction.getRightPredictionBorder(FILTER_REGION);
			
			if(border5 > border3){
				logger.warn("Safety border5 is bigger than safety border3!");
			}
			
			for (String me : prediction.getMobileMappings()) {
				if (knownMEs.containsKey(ref)
						&& knownMEs.get(ref).containsKey(me)) {
					for (String coordinates : knownMEs.get(ref).get(me)) {
						startKnown = Integer
								.parseInt(coordinates.split("-")[0]);
						endKnown = Integer.parseInt(coordinates.split("-")[1]);
						if ((border5 >= startKnown && border5 < endKnown)
								|| (border3 >= startKnown && border3 < endKnown)
								|| (border5 < startKnown && border3 >= endKnown)) {
							inKnownMe = true;
							if (me.equals("ALU")) {
								alu += 1;
							} else if (me.equals("SVA")) {
								sva += 1;
							} else if (me.equals("L1")) {
								l1 += 1;
							} else if (me.equals("HERV")) {
								herv += 1;
							}
							break;
						}
					}
					if(inKnownMe){
						break;
					}
				}
			}
			if(!inKnownMe){
				filteredPredictions.add(prediction);				
			}else{
				inKnownMe = false;
			}
		}
		logger.info("filtered alu: " + alu);
		logger.info("filtered sva: " + sva);
		logger.info("filtered l1: " + l1);
		logger.info("filtered herv (hervk): " + herv);
		return filteredPredictions;
	}

	
	/**
	 * This function reads in a RepeatMasker .out file and gets coordinates from the start
	 * or end positions (based on index param) of repetitive elements from this file and stores
	 * it in a HashMap.
	 * 
	 * @return Map containing chromosomes as keys and HashMap as values, containing
	 * mobile elements as keys and position start and position end as string values seperated by "-"
	 * @throws IOException
	 */
	public static Map<String, HashMap<String, HashSet<String>>> getKnownMEs() throws IOException{
		/*
		 * Map contains:
		 * - Chromosomes as key
		 * - Dictionary as value: with mobile element super family as key
		 *                        and position start (0-based inclusive) and position end (0-based exclusive)
		 *                        as sting value. Start and end are seperated by -
		 */
		Map<String, HashMap<String, HashSet<String>>> meMap = new HashMap<String, HashMap<String, HashSet<String>>>();
		InputStream is = null;
		String line = "";
		String chr = "";
		String element = "";
		String[] split;
		String startPos;
		String endPos;
		StringBuilder pos = new StringBuilder();
		boolean header = true;
		
		File f = new File(repmask_file);
		if (!f.exists() || !f.canRead())
			logger.equals(repmask_file + ": does not exist, or is not readable");
			
		is = new FileInputStream(f);
		
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		
		while ((line = br.readLine()) != null){
			if(!header){
				split = line.split("\t");
				chr = split[5];
				element = split[12].toUpperCase();
				startPos = split[6].trim();
				endPos = split[7].trim();
				pos.append(startPos);
				pos.append("-");
				pos.append(endPos);
				
				if(element.equals("OTHER")){
					element = "SVA";
				}else if(element.startsWith("ERV")){
					element = "HERV";
				}
				
				if (!meMap.containsKey(chr)){
					meMap.put(chr, new HashMap<String, HashSet<String>>());
					meMap.get(chr).put(element, new HashSet<String>());
					meMap.get(chr).get(element).add(pos.toString());
				}else if(!meMap.get(chr).containsKey(element)){
					meMap.get(chr).put(element, new HashSet<String>());
					meMap.get(chr).get(element).add(pos.toString());
				}else{
					meMap.get(chr).get(element).add(pos.toString());
				}
				pos.setLength(0);
			}
			else{
				header = false;
			}
		}
		br.close();
		is.close();

		
		return meMap;
	}
	
	
}


