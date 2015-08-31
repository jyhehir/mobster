package org.umcn.me.pairedend;

import net.sf.samtools.*;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.apache.commons.cli.*;
import org.umcn.gen.annotation.AnnotatedRegionSet;
import org.umcn.gen.region.LabeledRegion;
import org.umcn.gen.region.RegionComparator;
import org.umcn.me.sam.MateCluster;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.samexternal.SAMWriting;
import org.umcn.me.splitread.ClippedRead;
import org.umcn.me.splitread.InvalidHardClipCigarException;
import org.umcn.me.tabix.BlacklistAnnotation;
import org.umcn.me.tabix.RefGeneAnnotation;
import org.umcn.me.tabix.RepMaskAnnotation;
import org.umcn.me.tabix.SelfChainAnnotation;
import org.umcn.me.tabix.TabixBaseAnnotater;
import org.umcn.me.tabix.TabixReader;
import org.umcn.me.util.ClippedReadSet;
import org.umcn.me.util.CollectionUtil;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.ReadName;
import org.umcn.me.util.ReadNameOption;
import org.umcn.me.util.ReadNameRetriever;
import org.umcn.me.util.SampleBam;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.Vector;
import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;
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
	
	public static Logger logger = Logger.getLogger("AnchorClusterer");
	
	private static final int FILTER_REGION = 90;
	private static final String VERSION = "0.1.7b-GRIPS";
	
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
	
	private static Properties mobster_properties;
	
	
	public static String getVersionAndParameterInfo(String[] params){
		StringBuilder sb = new StringBuilder();
		Date date = new Date();
		sb.append("#Version: ");
		sb.append(VERSION);
		sb.append("\n");
		sb.append("#Command Line Initialization: ");
		
		for (String param : params){
			sb.append(param);
			sb.append(" ");
		}
		sb.append("\n");
		sb.append("#Creation date: ");
		sb.append(date.toString());
		sb.append("\n");

		return sb.toString();
	}
	
	public static String getVersionAndParameterInfo(Properties props){
		StringBuilder sb = new StringBuilder();
		Date date = new Date();
		sb.append("#Version: ");
		sb.append(VERSION);
		sb.append("\n");
		sb.append("#Properties file initialization, with following specified values : ");
		
		for (Object key : props.keySet()){
			String stringKey = (String) key;
			String value = props.getProperty(stringKey);
			sb.append(stringKey);
			sb.append("=");
			sb.append(value);
			sb.append(" ");
			
		}
		sb.append("\n");
		sb.append("#Creation date: ");
		sb.append(date.toString());
		sb.append("\n");
		
		return sb.toString();
	
	}
	
	
	public static void runFromPropertiesFile(Properties props) throws IOException{
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
		String commentHeader;
		String sample = ""; // default empty
		String outPrefix;
		int rpc; //default 2
		int overlap; // default 50 
		int maxdist; // default 450
		List<ReadName> anchorReads;
		List<ReadName> splitAnchorReads;
		Map<String, ReadName> readnameMap = new HashMap<String, ReadName>();
		
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
			repmask_file = props.getProperty(MobileDefinitions.REPEAT_MASK_FILE).trim();
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
		
		//TODO this is a copy and paste of the main function, needs refactoring
		
		logger.info("Using reads per cluster: " + rpc);
		logger.info("Using max overlap for double clusters: " + overlap);
		logger.info("Using max distance between 5 and 3 end clusters for double clusters of: " +
				maxdist);
		logger.info("using tmp: " + tmp);
		logger.info("using GRIPs mode: " + grips_mode);

		//Step 1: Clustering
		//Run simple clustering algorithm

		logger.info("Using simple clustering with search area of: " + search);
			
		runSimpleAnchorClusterer(anchor, anchorIndex, clusterBam, search, rpc);

		
		//Step 3: Make single & double cluster predictions out of the clusters
		Vector<MobilePrediction> matePredictions;
		
		try {
			matePredictions = clusterMateClusters(clusterBam, new File(clusterBam.toString().replaceAll(".bam$", ".bai")),
					overlap, maxdist);
			
			//Step 4: do split clustering when available
			if(splitAnchor != null){
				logger.info("Number of initial minimum hits for split clusters:" + min_initial_cluster_size_split);
				clusterSplitAnchorsToBAM(splitAnchor, splitAnchorIndex,
						new File(outPrefix + "_splitcluster.bam"), min_initial_cluster_size_split);
				
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
			
			if (mobster_properties.containsKey(MobileDefinitions.FILTER_OVERLAPPING_PREDICTIONS) &&
					Boolean.parseBoolean(mobster_properties.getProperty(MobileDefinitions.FILTER_OVERLAPPING_PREDICTIONS))){
				logger.info("Anchorclusterer: will remove overlapping predictions");
				matePredictions = GRIPFunctions.removeOverlappingPredictions(matePredictions);
			}
			commentHeader = getVersionAndParameterInfo(props);
			if (grips_mode){
				writeGRIPSToFile(outPrefix + "_GRIPS_unfiltered.txt", matePredictions, commentHeader, false);
			}
			
			matePredictions = filterByMinTotalHits(matePredictions, min_total_hits);
			
			if (!grips_mode){
				matePredictions = filterKnownMEs(getKnownMEs(), matePredictions);
			}
			
			//Filter for multiple occuring source genes if detection is done using GRIPS
			if (props.containsKey(MobileDefinitions.GRIPS_MAX_SOURCE_GENES)){
				int maxSourceGenes = Integer.parseInt(props.getProperty(MobileDefinitions.GRIPS_MAX_SOURCE_GENES));
				logger.info("Filtering GRIPS for max source genes: " + maxSourceGenes);
				matePredictions = GRIPFunctions.reducePredictionsBasedOnSource(matePredictions, maxSourceGenes);
			}
			
			
			if (grips_mode){
				//Prestep 1 for GRIPS: extracting the read names
				anchorReads = extractReadnames(anchor);
				splitAnchorReads = extractReadnames(splitAnchor);
				
				
				System.out.println("Anchor reads size: " + anchorReads.size());
				System.out.println("Split anchor reads size: " + splitAnchorReads.size());
				readnameMap = CollectionUtil.readNamesToMap(anchorReads);
				System.out.println("read name map size 1: " + readnameMap.size());
				readnameMap.putAll(CollectionUtil.readNamesToMap(splitAnchorReads));
				System.out.println("read name map size 2: " + readnameMap.size());
				
				
				//Now only keep the reads in map which made it to the clusters supporting the predictions
				Set<String> reads = new HashSet<String>();
				for (MobilePrediction pred : matePredictions){
					reads.addAll(pred.getReadNamesFromMateClusters());
					reads.addAll(pred.getReadNamesFromSplitClusters());
				}
				System.out.println("Number of total supporting reads from clusters: " + reads.size());
				Set<String> keys = readnameMap.keySet();
				keys.retainAll(reads);
				System.out.println("Number of keys in map: " + readnameMap.size());
				
				//annotate the predictions
				String refGeneLoc = props.getProperty(MobileDefinitions.GRIPS_TRANSCRIPT);
				String repMaskLoc = props.getProperty(MobileDefinitions.GRIPS_REPMASK);
				String blackListLoc = props.getProperty(MobileDefinitions.GRIPS_BLACKLIST);
				String selfChainLoc = props.getProperty(MobileDefinitions.GRIPS_SELFCHAIN);
				
				try {
					annotateRefGene(readnameMap.values(), refGeneLoc);
					annotateRepMask(readnameMap.values(), repMaskLoc, true);
					annotateBlacklist(readnameMap.values(), blackListLoc);
					annotateSelfChain(readnameMap.values(), selfChainLoc);
				} catch (java.text.ParseException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				
				//Add in the repmask / refgene annotation
				for (MobilePrediction pred : matePredictions){
					pred.parseReadNames(readnameMap);
				}
				writeGRIPSToFile(outPrefix + "_GRIPS_predictons.txt", matePredictions, commentHeader, false);
				writeGRIPSToFile(outPrefix + "_GRIPS_MORECONFIDENT_predictions.txt", matePredictions, commentHeader, true);
			}else{
				writePredictionsToFile(outPrefix + "_predictions.txt", matePredictions, commentHeader, sample);
			}

			long end = System.currentTimeMillis();
			long millis = end - start;
			String time = String.format("%d min, %d sec", 
				    TimeUnit.MILLISECONDS.toMinutes(millis),
				    TimeUnit.MILLISECONDS.toSeconds(millis) - 
				    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
				);
			
			logger.info("Anchorclusterer ran in : " + time);
		} catch (IOException e) {
			logger.error("Anchorclusterer: Error in opening files: " + e.getMessage());
		} finally {
			if (tmp != null && ! tmp.delete() ){
				logger.error("Anchor Clusterer: Could not delete temp: " + tmp);
			}
		}

		
	}

	public static void annotateSelfChain(Collection<ReadName> reads, String chainLocation) throws IOException, java.text.ParseException{
		TabixBaseAnnotater<SelfChainAnnotation> tba = new TabixBaseAnnotater<SelfChainAnnotation>(chainLocation, new SelfChainAnnotation());
		
		for (ReadName name : reads){
			if ( name.mateIsMapped ){
				name.setSelfChain(tba.queryOverlapping(name.toPositionString())); //intentionally query from the anchor position
			}
			
		}
		
	}
	
	public static void annotateRefGene(Collection<ReadName> reads, String transcriptLocation)
			throws IOException, java.text.ParseException {
		TabixBaseAnnotater<RefGeneAnnotation> tba = new TabixBaseAnnotater<RefGeneAnnotation>(transcriptLocation, new RefGeneAnnotation());
		
		
		for (ReadName name : reads){
			if( name.mateIsMapped){
				name.setMateRefGeneAnnotation(tba.queryOverlapping(name.mateToPositionString()));
			}
			if (name.isMapped){
				name.setRefGeneAnnotation(tba.queryOverlapping(name.toPositionString()));
			}
		}
	}
	
	public static void annotateBlacklist(Collection<ReadName> reads, String blacklistLocation) throws IOException, java.text.ParseException{
		TabixBaseAnnotater<BlacklistAnnotation> tba = new TabixBaseAnnotater<BlacklistAnnotation>(blacklistLocation, new BlacklistAnnotation());
		
		for (ReadName name : reads){
			if (name.isMapped){
				name.setBlacklist(tba.queryOverlapping(name.toPositionString()));
			}
			if (name.mateIsMapped){
				name.setMateBlacklist(tba.queryOverlapping(name.mateToPositionString()));
			}
		}
		
	}
	
	public static void annotateRepMask(Collection<ReadName> reads, String repmaskLocation, boolean annotateMates)
			throws IOException, java.text.ParseException {
		TabixBaseAnnotater<RepMaskAnnotation> tba = new TabixBaseAnnotater<RepMaskAnnotation>(repmaskLocation, new RepMaskAnnotation());
		
		for (ReadName name : reads){
			if(name.isMapped){
				name.setRepMaskAnnotation(tba.queryOverlapping(name.toPositionString()));
			}
			if(annotateMates && name.mateIsMapped){
				name.setMateRepMaskAnnotation(tba.queryOverlapping(name.mateToPositionString()));
			}
		}
	}

	public static List<ReadName> extractReadnames(File anchor) {
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).autoPrefixReference(true).build();
		ReadNameRetriever retriever = new ReadNameRetriever(anchor, option);
		List<ReadName> reads = new ArrayList<ReadName>();
		
		for (ReadName read : retriever){
			reads.add(read);
		}
		return reads;
	}
	
	public static void main(String[] args) {
		Options options;
		File anchor;
		File anchorIndex;
		File clusterBam;
		File splitAnchor = null;
		File splitAnchorIndex = null;
		String commentHeader;
		String sample;
		int rpc;
		int overlap;
		int maxdist;
		Properties prop = new Properties();
		
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
				
				if (line.hasOption("properties")){
					prop.load(new FileInputStream(line.getOptionValue("properties")));
					runFromPropertiesFile(prop);
					System.exit(0);
				}
				
				anchor = new File(line.getOptionValue("in"));
				anchorIndex = new File(line.getOptionValue("in").replaceAll(".bam$", ".bai"));
				clusterBam = new File(line.getOptionValue("out") + "_clusters.bam");
				sd_frag_size = Integer.parseInt(line.getOptionValue("sd"));
				mean_frag_size = Integer.parseInt(line.getOptionValue("mfl"));
				rpc = Integer.parseInt(line.getOptionValue("rpc", Integer.toString(rpc_default)));
				search = Integer.parseInt(line.getOptionValue("search", Integer.toString(search)));
				overlap = Integer.parseInt(line.getOptionValue("overlap", Integer.toString(overlap_default)));
				maxdist = Integer.parseInt(line.getOptionValue("maxdist", Integer.toString(maxdist_default)));
				min_total_hits = Integer.parseInt(line.getOptionValue("mintotal", Integer.toString(min_total_hits)));
				min_initial_cluster_size_split = Integer.parseInt(line.getOptionValue("splithits", Integer.toString(min_initial_cluster_size_split)));
				//samdir = line.getOptionValue("samdir", "");
				//lenient_search = line.hasOption("lenient");
				sample = line.getOptionValue("sample", "");
				multiple_sample_calling = line.hasOption("multiplesample");
				filter_by_read_counts_single_sample = line.hasOption("multisample_stringent");
				percentile_99_fragment = Integer.parseInt(line.getOptionValue("maxclust", Integer.toString(percentile_99_fragment)));
				TMP = line.getOptionValue("tmp", System.getProperty("java.io.tmpdir"));
				memory = Integer.parseInt(line.getOptionValue("max_memory", Integer.toString(SAMWriting.MAX_RECORDS_IN_RAM)));				
				grips_mode = line.hasOption("grips");
				
				if(line.hasOption("insplit")){
					splitAnchor = new File(line.getOptionValue("insplit"));
					splitAnchorIndex = new File(line.getOptionValue("insplit").replaceAll(".bam$", ".bai"));
				}
				
				repmask_file = line.getOptionValue("repmask", repmask_file);
				

				logger.info("Using reads per cluster: " + rpc);
				logger.info("Using max overlap for double clusters: " + overlap);
				logger.info("Using max distance between 5 and 3 end clusters for double clusters of: " +
						maxdist);
				logger.info("Using GRIPs mode: " + grips_mode);
				
				//Step 1: Clustering
				//Run simple clustering algorithm

				logger.info("Using simple clustering with search area of: " + search);
					
				runSimpleAnchorClusterer(anchor, anchorIndex, clusterBam, search, rpc);

				
				//Step 2: coordinate sorting the clusters (each cluster is a
				//read name in the bam file)
//				execCommand(samdir + "samtools sort " + clusterBam.toString() + " " +
//						line.getOptionValue("out") + "_clusters_sorted");
//				
//				execCommand(samdir + "samtools index " + line.getOptionValue("out") + "_clusters_sorted.bam");
				
				//Step 3: Make single & double cluster predictions out of the clusters
				Vector<MobilePrediction> matePredictions;
				
				matePredictions = clusterMateClusters(clusterBam, new File(clusterBam.toString().replaceAll(".bam$", ".bai")),
						overlap, maxdist);
				
				//Step 4: do split clustering when available
				if(splitAnchor != null){
					logger.info("Number of initial minimum hits for split clusters:" + min_initial_cluster_size_split);
					clusterSplitAnchorsToBAM(splitAnchor, splitAnchorIndex,
							new File(line.getOptionValue("out") + "_splitcluster.bam"), min_initial_cluster_size_split);
					
//					execCommand(samdir + "samtools sort " + line.getOptionValue("out") + "_splitcluster.bam" +
//					 " " + line.getOptionValue("out") + "_splitcluster_sorted");
//					execCommand(samdir + "samtools index " + line.getOptionValue("out") + "_splitcluster_sorted.bam");
					
					try {
						matePredictions = mergeSplitAndMateClusters(matePredictions, new File(line.getOptionValue("out") +
								"_splitcluster.bam"), new File(line.getOptionValue("out") + "_splitcluster.bai"));
					} catch (IllegalSAMPairException e) {
						logger.error(e.getMessage());
					}
				}
				
				matePredictions = mergePredictions(matePredictions, overlap, maxdist);
				
				matePredictions = filterByMinTotalHits(matePredictions, min_total_hits);
				matePredictions = filterKnownMEs(getKnownMEs(), matePredictions);
				commentHeader = getVersionAndParameterInfo(args);
				writePredictionsToFile(line.getOptionValue("out") + "_predictions.txt", matePredictions, commentHeader, sample);
				
				long end = System.currentTimeMillis();
				long millis = end - start;
				String time = String.format("%d min, %d sec", 
					    TimeUnit.MILLISECONDS.toMinutes(millis),
					    TimeUnit.MILLISECONDS.toSeconds(millis) - 
					    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
					);
				
				logger.info("Anchorclusterer ran in : " + time);
				
				
			}  catch (IOException e){
				logger.error("Error in opening files: " + e.getMessage());
			} catch (org.apache.commons.cli.ParseException e) {
				// TODO Auto-generated catch block
				logger.error("Error in parsing CLI arguments: " + e.getMessage());
			}
			
		}

	}
	
	private static Vector<MobilePrediction> filterByMinTotalHits(
			Vector<MobilePrediction> matePredictions, int min_total_hits2) {
		
		Vector<MobilePrediction> copy_preds = new Vector<MobilePrediction>();
		int c = 0;
		boolean pass;
				
		for (MobilePrediction pred : matePredictions){
			
			pass = false;
			
			if (!multiple_sample_calling || !filter_by_read_counts_single_sample){
				if (pred.getLeftTotalHits() + pred.getRightTotalHits() >= min_total_hits2){
					copy_preds.add(pred);
					pass = true;
				}
			}
			
			else{
				Map<String, Integer> sampleCounts = pred.getSampleCounts();
				for (String sample : sampleCounts.keySet()){
					if (sampleCounts.get(sample) >= min_total_hits2){
						copy_preds.add(pred);
						pass = true;
						break;
					}
				}
			}
			if (!pass){
				c++;
			}
		}
		
		logger.info(c +" predictions did not meet the minimum required read support of: " + min_total_hits2 );
		return copy_preds;
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
		Map<LabeledRegion, MobilePrediction> regionMap = CollectionUtil.mobilePredictionsToRegions(predictions, overlap_default);
		AnnotatedRegionSet<LabeledRegion> regions = new AnnotatedRegionSet<LabeledRegion>(regionMap.keySet());
		Set<MobilePrediction> uniquePredictions = new HashSet<MobilePrediction>();
		
		regions.sort();
		int c = 0;
		RegionComparator<LabeledRegion, LabeledRegion> rc = new RegionComparator<LabeledRegion, LabeledRegion>();
		for (LabeledRegion r : regions){
			AnnotatedRegionSet<LabeledRegion> overlap = rc.getOverlapRegionsFast(r, regions);
			for (LabeledRegion or : overlap){
				if (or.equals(r)){
					continue;
				}else{
					if ( mergeMobilePredictions(r, or, regionMap, overlapdef, maxdist)){
						c++;
					}
				}
			}
		}
		
		for (LabeledRegion r : regions){
			uniquePredictions.add(regionMap.get(r));
		}
		
		logger.info("Number of merged Mobile Predictions: " + c);
		logger.info("Number of predictions before merging: " + predictions.size());
		logger.info("Number of predictions after merging: " + uniquePredictions.size());
		
		
		Map<LabeledRegion, MobilePrediction> regionMap2 = CollectionUtil.mobilePredictionsToRegions(new Vector<MobilePrediction>(uniquePredictions), overlap_default);
		AnnotatedRegionSet<LabeledRegion> regions2 = new AnnotatedRegionSet<LabeledRegion>(regionMap2.keySet());
		regions2.sort();
		Vector<MobilePrediction> mergedPredictions =  new Vector<MobilePrediction>();
		for (LabeledRegion r : regions2){
			mergedPredictions.add(regionMap2.get(r));
		}
		
		return mergedPredictions;
		
	}
	
	private static boolean mergeMobilePredictions(LabeledRegion r1, LabeledRegion r2,
			Map<LabeledRegion, MobilePrediction> labeledToMobile, int overlapdef, int maxdist){
		
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
	
	public static void writePredictionsToFile(String outString, Vector<MobilePrediction> predictions, String comment, String sampleName){
		try {
			PrintWriter outFile = new PrintWriter(new FileWriter(outString), true);
			outFile.print(comment);
			
			outFile.print(predictions.get(0).getHeader() + "\n");
			
			for(MobilePrediction pred : predictions){
				//pred.setSampleName(sampleName);
				outFile.print(pred + "\n");
			}
			
			outFile.close();
		} catch (IOException e) {
			logger.error("Failed to write prediction file: " + e.getMessage());
		}
		
		
	}
	
	public static void writeGRIPSToFile(String outString, Vector<MobilePrediction> predictions, String comment, boolean filter){
		try {
			PrintWriter outFile = new PrintWriter(new FileWriter(outString), true);
			
			boolean writtenHeader = false;
			outFile.print(comment);
			for(MobilePrediction pred : predictions){
				if (! writtenHeader){
					outFile.println(pred.toGripsHeader());
					writtenHeader = true;
				}

				if (filter && ! pred.gripNeedsFiltering()){
					outFile.println(pred.toGripsString());
				}else if (! filter){
					outFile.println(pred.toGripsString());
				}

			}
			
			outFile.close();
		} catch (IOException e) {
			logger.error("Failed to write prediction file: " + e.getMessage());
		}
		
		
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
	public static void runSimpleAnchorClusterer(File anchor, File index, File clusterBam, int searchRegion, int minReads){

		SAMFileReader input = new SAMSilentReader(anchor, index);
		SAMFileHeader header = input.getFileHeader();
		SampleBam sampleCalling = SampleBam.SINGLESAMPLE;
		
		header = editSAMFileHeader(header);
		SAMFileWriter outputSam = SAMWriting.makeSAMWriter(clusterBam, header, new File(TMP), memory, SAMFileHeader.SortOrder.coordinate);
		
		Vector<MateCluster<SAMRecord>> mobileClusters = new Vector<MateCluster<SAMRecord>>();
		
		int c = 0;
		int skippedClustersBecauseOfNotSameRefMapping = 0;
		double minPercentSameMateRefMapping = 0.0;
		int maxDiffMateMapping = Integer.MAX_VALUE;
		
		if (mobster_properties.containsKey(MobileDefinitions.GRIPS_MIN_PERCENT_SAME_REF_MATE_MAPPING)){
			minPercentSameMateRefMapping = Double.parseDouble(mobster_properties.getProperty(MobileDefinitions.GRIPS_MIN_PERCENT_SAME_REF_MATE_MAPPING));
		}
		
		if (mobster_properties.containsKey(MobileDefinitions.GRIPS_MAX_DIFF_MATE_MAPPINGS)){
			maxDiffMateMapping = Integer.parseInt(mobster_properties.getProperty(MobileDefinitions.GRIPS_MAX_DIFF_MATE_MAPPINGS));
		}
		
		logger.info("Using a minPercentSameMateRefMapping threshold of: " + minPercentSameMateRefMapping);
		logger.info("Using a max different mate chr threshold of: " + maxDiffMateMapping);
		
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
						currentCluster.writeClusterToSAMWriter(outputSam, Integer.toString(c));
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
				cluster.writeClusterToSAMWriter(outputSam, Integer.toString(c));
			}
		}
		logger.info("Number of skipped clusters because mates do not align to same chromosome:" + skippedClustersBecauseOfNotSameRefMapping);
		logger.info("Found nr of clusters: " + c);

		input.close();
		outputSam.close();
	}
	
	public static void clusterSplitAnchorsToBAM(File anchor, File index, File outBam, int minClustersize){
		
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
		header = editSAMFileHeader(header);
		header2 = editSAMFileHeader(header2);
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
	
	
	
	private static SAMFileHeader editSAMFileHeader(SAMFileHeader header) {
		List<SAMSequenceRecord> sequenceList = new ArrayList<SAMSequenceRecord>();
		
		for ( SAMSequenceRecord sequenceRec : header.getSequenceDictionary().getSequences()){
			StringBuilder sequenceBuilder = new StringBuilder();
			SAMSequenceRecord newSeq;
			if(!sequenceRec.getSequenceName().startsWith("chr")){
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
	 * @param predictionFile name of output prediction .txt file
	 * @throws IOException
	 */
	public static Vector<MobilePrediction> clusterMateClusters(File clusterIn, File clusterIndex, int overlap, int maxdist) throws IOException{
		
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
	 * Write a tab delimited file
	 * 
	 * @param lines Nested vector. Each inner vector contains a string value for each column.
	 * @param header Vector. Each String inside vector contains a column header.
	 * @param outname Name of outputted file.
	 * @param delim Delimiter to be used between column names and values.
	 * @throws IOException
	 */
	public static void writeDelimitedFile(Vector<Vector<String>> lines, Vector<String> header, String outname,
			String delim)
		throws IOException{
		
		PrintWriter outFile = new PrintWriter(new FileWriter(outname), true);
		writeDelimitedLine(header, outFile, delim);
		
		for(Vector<String> line : lines){
			writeDelimitedLine(line, outFile, delim);
		}
		
		outFile.close();
		
		
	}
	
	/**
	 * Write a single delimited line from a PrintWriter object.
	 * @param line Vector containing string values to be printed.
	 * @param out PrintWriter object to write from.
	 * @param delim Delimiter used for seperating string values in parameter line.
	 */
	public static void writeDelimitedLine(Vector<String> line, PrintWriter out, String delim){
		
		for(int i = 0; i < line.size(); i++){
			out.print(line.get(i));
			if (i != line.size() - 1){
				out.print(delim);
			}
		}
		out.println();
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
		
		is = new FileInputStream(repmask_file);
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
