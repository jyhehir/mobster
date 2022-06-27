package org.umcn.me.pairedend;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.tabix.BlacklistAnnotation;
import org.umcn.me.tabix.RefGeneAnnotation;
import org.umcn.me.tabix.RepMaskAnnotation;
import org.umcn.me.util.*;

import java.util.*;

public class MobilePrediction {

	protected static final int MIN_DISCORDANT_UX = 2;

	protected String original_reference = "";
	protected ReferenceGenome reference_genome;
	protected boolean output_somatic = false;

	protected int median_fragment_length;
	protected int sd_fragment_length;

	protected int left_mate_cluster_border = 0;
	protected int left_cluster_length = 0;
	protected int left_mate_hits = 0;
	protected int left_aligned_split_hits = 0; //aligned on left thus clipped on right
	protected int left_aligned_split_border = 0;
	protected int left_aligned_polyA_hits = 0;
	protected int left_aligned_polyT_hits = 0;

	protected int right_mate_cluster_border = 0;
	protected int right_cluster_length = 0;
	protected int right_mate_hits = 0;
	protected int right_aligned_split_hits = 0; //aligned on right thus clipped on left
	protected int right_aligned_split_border = 0;
	protected int right_aligned_polyA_hits = 0;
	protected int right_aligned_polyT_hits = 0;

	protected int unique_hits = 0;
	protected int multiple_hits = 0;
	protected int unmapped_hits = 0;

	protected int max_expected_cluster_size;


	protected String leftclipped_max_distance = "-1";
	protected String rightclipped_max_distance = "-1";
	protected String leftclipped_same_fraction = "-1";
	protected String rightclipped_same_fraction = "-1";
	protected String clipped_avg_qual = "-1";
	protected String clipped_avg_len = "-1";

	protected boolean merged = false;

	protected String mate_mobile_mapping = "";
	protected Set<String> mobile_mappings = new HashSet<String>();
	protected Set<String> sample_names = new HashSet<String>();
	protected Set<String> split_read_names = new HashSet<String>();
	protected Set<String> discordant_read_names = new HashSet<String>();
	protected Map<String, Integer> sample_counts = new HashMap<String, Integer>();
	protected Map<String, Integer> non_supporting_sample_counts = new HashMap<String, Integer>();;
	protected Map<String, Double> sample_vafs = new HashMap<String, Double>();;
	protected String somatic = "";
	
	protected Map<String, Integer> refseqMateCounts = new HashMap<String, Integer>();
	protected Map<String, Integer> repMaskAnchor_counts = new HashMap<String, Integer>();
	protected Map<String, Integer> repMaskMate_counts = new HashMap<String, Integer>();
	protected Set<String> refSeqAnchorLocations = new HashSet<String>();
	protected Map<String, Integer> blacklistAnchorCounts = new HashMap<String, Integer>();
	protected Map<String, Integer> blacklistMateCounts = new HashMap<String, Integer>();
	protected Map<String, Integer> repFamilyAnchorCounts = new HashMap<String, Integer>();
	protected Map<String, Integer> repFamilyMateCounts = new HashMap<String, Integer>();
	protected List<Double> chainScores = new ArrayList<Double>();
	protected List<Integer> insertSizes = new ArrayList<Integer>();

	protected SAMRecord samrecord_splitcluster = null;

	protected SAMRecord right_cluster_sam_record;

	protected SAMRecord left_cluster_sam_record;

	protected int selfChainOverlap = 0;
	protected String selfChainOverlapString = "";

	protected Map<String, Integer> refSeqMateUUCounts = new HashMap<String, Integer>();

	protected int left_mate_median_mapq = -1;

	protected int right_mate_median_mapq = -1;


	public static Logger logger = Logger.getLogger(AnchorClusterer.class.getName());

	public final static int SINGLE_CLUSTER_BORDER_FUZZINESS = 20;
	public final static int EXPECTED_DUPLICATION_LENGTH = 15;

	public final static String DUPLICATION = "duplication";
	public final static String DELETION = "deletion";
	public final static String NO_TSD = "no_tsd";
	public final static String UNKNOWN = "unknown";

	public final static String SOMATIC = "somatic";
	public final static String GERMLINE = "germline";
	public final static String ARTIFACT = "artifact";

	public final static String COLUMN_REFERENCE = "Chr";
	public final static String COLUMN_MOBILE = "Mobile Element";
	public final static String COLUMN_INSERTPOINT = "Insert Point";
	public final static String COLUMN_BORDER5 = "border5";
	public final static String COLUMN_BORDER3 = "border3";
	public final static String COLUMN_ENDINSERTPOINT = "End Insert Point";
	public final static String COLUMN_ENDBORDER5 = "end border5";
	public final static String COLUMN_ENDBORDER3 = "end border3";
	public final static String COLUMN_MERGED = "merged";
	public final static String COLUMN_CLUSTER5_LEN = "cluster5 length";
	public final static String COLUMN_CLUSTER3_LEN = "cluster3 length";
	public final static String COLUMN_CLUSTER5_HITS = "cluster5 hits";
	public final static String COLUMN_CLUSTER3_HITS = "cluster3 hits";
	public final static String COLUMN_SPLIT5_HITS = "split5 hits";
	public final static String COLUMN_SPLIT3_HITS = "split3 hits";
	public final static String COLUMN_POLYA5_HITS = "polyA5 hits";
	public final static String COLUMN_POLYT5_HITS = "polyT5 hits";
	public final static String COLUMN_POLYA3_HITS = "polyA3 hits";
	public final static String COLUMN_POLYT3_HITS = "polyT3 hits";

	public final static String COLUMN_UNIQUE_HITS = "original discordant unique";
	public final static String COLUMN_MULTIPLE_HITS = "original multiple";
	public final static String COLUMN_UNMAPPED_HITS = "original unmapped";
	public final static String COLUMN_LEFTCLIPPED_MAX_DISTANCE = "leftclipped max dist";
	public final static String COLUMN_RIGHTCLIPPED_MAX_DISTANCE = "rightclipped max dist";
	public final static String COLUMN_LEFTCLIPPED_FRAC_DISTANCE = "leftclipped same pos";
	public final static String COLUMN_RIGHTCLIPPED_FRAC_DISTANCE = "rightclipped same pos";
	public final static String COLUMN_CLIPPED_AVG_QUAL = "clipped avg qual";
	public final static String COLUMN_AVG_CLIPPED_LEN = "clipped avg length";
	public final static String COLUMN_TSD = "target site duplication";
	public final static String COLUMN_TSD_LEN = "target site duplication length";
	public final static String COLUMN_TSD_SEQ = "target site duplication sequence";
	public final static String COLUMN_SAMPLE = "sample";
	public final static String COLUMN_SAMPLE_COUNT = "sample_counts";
	public final static String COLUMN_NON_SUPPORTING_SAMPLE_COUNT = "non_supporting_sample_counts";
	protected static final String COLUMN_VAF = "variant allele fraction";
	protected static final String COLUMN_SOMATIC = "somatic";

	protected final static String[] HEADER = {COLUMN_REFERENCE, COLUMN_MOBILE, COLUMN_INSERTPOINT, COLUMN_BORDER5, COLUMN_BORDER3, COLUMN_ENDINSERTPOINT,
											COLUMN_ENDBORDER5, COLUMN_ENDBORDER3, COLUMN_MERGED, COLUMN_SAMPLE, COLUMN_SAMPLE_COUNT, COLUMN_NON_SUPPORTING_SAMPLE_COUNT,
											COLUMN_CLUSTER5_LEN, COLUMN_CLUSTER3_LEN, COLUMN_CLUSTER5_HITS, COLUMN_CLUSTER3_HITS, COLUMN_SPLIT5_HITS,
											COLUMN_SPLIT3_HITS, COLUMN_POLYA5_HITS, COLUMN_POLYT5_HITS, COLUMN_POLYA3_HITS, COLUMN_POLYT3_HITS,
											COLUMN_UNIQUE_HITS, COLUMN_MULTIPLE_HITS, COLUMN_UNMAPPED_HITS, COLUMN_LEFTCLIPPED_MAX_DISTANCE, COLUMN_RIGHTCLIPPED_MAX_DISTANCE,
	                                        COLUMN_LEFTCLIPPED_FRAC_DISTANCE, COLUMN_RIGHTCLIPPED_FRAC_DISTANCE, COLUMN_CLIPPED_AVG_QUAL, COLUMN_AVG_CLIPPED_LEN, COLUMN_TSD,
											COLUMN_TSD_LEN, COLUMN_TSD_SEQ, COLUMN_VAF, COLUMN_SOMATIC};

	protected static final int MIN_INSERT_SIZE = 10000;

//	protected final static String[] GRIPS_HEADER = {COLUMN_GRIP_REFSEQ, COLUMN_GRIP_REPNAME_ANCHOR, COLUMN_GRIP_REPCLASS_ANCHOR,
//												  COLUMN_GRIP_REPFAMILY_ANCHOR, COLUMN_GRIP_REPNAME_MATE, COLUMN_GRIP_REPCLASS_MATE,
//												  COLUMN_GRIP_REPFAMILY_MATE};

	protected static Map<String, String> features = new HashMap<String, String>();

	public MobilePrediction(int mfl, int sd, int maxExpectedClusterSize){
		this.median_fragment_length = mfl;
		this.sd_fragment_length = sd;
		this.max_expected_cluster_size = maxExpectedClusterSize;
	}

	public MobilePrediction(int mfl, int sd, int maxExpectedClusterSize, SAMRecord cluster) throws IllegalSAMPairException{
		this(mfl, sd, maxExpectedClusterSize);
		parseSAMRecordCluster(cluster);
	}

	public MobilePrediction(int mfl, int sd, int maxExpectedClusterSize, SAMRecord cluster1, SAMRecord cluster2)
			throws IllegalSAMPairException{
		this(mfl, sd, maxExpectedClusterSize, cluster1);
		parseSAMRecordCluster(cluster2);
	}

	public boolean addDiscordantReadsOfMobilePrediction(MobilePrediction pred, int overlap, int maxdist) throws IllegalSAMPairException{

		boolean success = false;

		overlap = - Math.abs(overlap);//make sure overlap is negative


		//dont merge mobile predictions of different families
		if (!pred.mobile_mappings.equals(this.mobile_mappings)){
			return false;
		}

		if (!pred.original_reference.equals(this.original_reference)){
			return false;
		}

		SAMRecord leftClusterToAdd = pred.left_cluster_sam_record;
		SAMRecord rightClusterToAdd = pred.right_cluster_sam_record;

		if (pred.onlyLeftDiscordantReads() && this.onlyRightDiscordantReads()){
			if (this.right_mate_cluster_border - pred.left_mate_cluster_border >= maxdist){
				return false;
			}
		}

		if (pred.onlyRightDiscordantReads() && this.onlyLeftDiscordantReads()){
			if (pred.right_mate_cluster_border - this.left_mate_cluster_border >= maxdist){
				return false;
			}
		}

		if (leftClusterToAdd != null &&
				!this.discordant_read_names.contains(leftClusterToAdd.getReadName())){

			if (this.right_cluster_sam_record == null || (this.right_cluster_sam_record != null && this.right_cluster_sam_record.getAlignmentStart() -
					leftClusterToAdd.getAlignmentEnd() >= overlap)){


				if (this.left_cluster_sam_record != null){
					int[] coordinates = {this.left_cluster_sam_record.getAlignmentStart(), this.left_cluster_sam_record.getAlignmentEnd(),
							leftClusterToAdd.getAlignmentStart(), leftClusterToAdd.getAlignmentEnd()};
					Arrays.sort(coordinates);

					if (coordinates[2] - coordinates[1] >= maxdist){
						logger.info("Same orientation clusters are too distant to be merged in region: " + this.getChromosome() + ":" +
								this.getLeftPredictionBorder() + "-" + this.getEndRightPredictionBorder());
						return false;
					}

					this.left_cluster_length = coordinates[3] - coordinates[0] + 1;
					this.left_mate_cluster_border = coordinates[3];
				}else{
					this.left_cluster_sam_record = leftClusterToAdd;
					this.left_cluster_length = leftClusterToAdd.getAlignmentEnd() - leftClusterToAdd.getAlignmentStart() + 1;
					this.left_mate_cluster_border = leftClusterToAdd.getAlignmentEnd();
				}

				this.unique_hits += Integer.parseInt(leftClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_UNIQUE_HITS).toString());
				this.multiple_hits += Integer.parseInt(leftClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_MULTIPLE_HITS).toString());
				this.unmapped_hits += Integer.parseInt(leftClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_UNMAPPED_HITS).toString());
				this.discordant_read_names.add(leftClusterToAdd.getReadName());
				this.left_mate_hits += Integer.parseInt(leftClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_CLUSTER_HITS).toString());
				this.parseSampleCountTag(leftClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_SAMPLECOUNT).toString());

				this.merged = true;
				success = true;

			} else{
				logger.info("Did not merge in following region because left cluster to add had more than allowed " +
						"overlap with this.right_mate_cluster in region: " + this.getChromosome() + ":" +
						this.getLeftPredictionBorder() + "-" + this.getEndRightPredictionBorder());
				return false;
			}


		}

		if(rightClusterToAdd != null &&
				!this.discordant_read_names.contains(rightClusterToAdd.getReadName())){
			if (this.left_cluster_sam_record == null || (this.left_cluster_sam_record != null && rightClusterToAdd.getAlignmentStart() -
					this.left_cluster_sam_record.getAlignmentEnd() >= overlap)){


				if (this.right_cluster_sam_record != null){
					int[] coordinates = {this.right_cluster_sam_record.getAlignmentStart(), this.right_cluster_sam_record.getAlignmentEnd(),
						rightClusterToAdd.getAlignmentStart(), rightClusterToAdd.getAlignmentEnd()};
					Arrays.sort(coordinates);
					if (coordinates[2] - coordinates[1] >= maxdist){
						logger.info("Same orientation clusters are too distant to be merged in region: " + this.getChromosome() + ":" +
								this.getLeftPredictionBorder() + "-" + this.getEndRightPredictionBorder());
						return false;
					}
					this.right_cluster_length = coordinates[3] - coordinates[0] + 1;
					this.right_mate_cluster_border = coordinates[0];
				}else{
					this.right_cluster_sam_record = rightClusterToAdd;
					this.right_cluster_length = rightClusterToAdd.getAlignmentEnd() - rightClusterToAdd.getAlignmentStart() + 1;
					this.right_mate_cluster_border = rightClusterToAdd.getAlignmentStart();
				}

				this.unique_hits += Integer.parseInt(rightClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_UNIQUE_HITS).toString());
				this.multiple_hits += Integer.parseInt(rightClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_MULTIPLE_HITS).toString());
				this.unmapped_hits += Integer.parseInt(rightClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_UNMAPPED_HITS).toString());
				this.discordant_read_names.add(rightClusterToAdd.getReadName());
				this.right_mate_hits += Integer.parseInt(rightClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_CLUSTER_HITS).toString());
				this.parseSampleCountTag(rightClusterToAdd.getAttribute(MobileDefinitions.SAM_TAG_SAMPLECOUNT).toString());
				success = true;
				this.merged = true;
			}else{
				logger.info("Did not merge in following region because right cluster to add had more than allowed " +
						"overlap with this.left_mate_cluster in region: " + this.getChromosome() + ":" +
						this.getLeftPredictionBorder() + "-" + this.getEndRightPredictionBorder());
				return false;
			}
		}

		//Add extra split cluster when this prediction does not have any split clusters
		if(pred.samrecord_splitcluster != null && this.samrecord_splitcluster == null){
			logger.info("Added split cluster when merging in region: " + this.getChromosome() + ":" +
					this.getLeftPredictionBorder() +  "-" + this.getEndRightPredictionBorder());
			this.parseSAMRecordCluster(pred.samrecord_splitcluster);
			success = true;
			this.merged = true;
		}

		return success;

	}

	public void parseReadNames(Map<String,ReadName> readnameMap){

		List<String> mateClusterReadNames = this.getReadNamesFromMateClusters();
		List<String> splitClusterReadNames = this.getReadNamesFromSplitClusters();

		//Step 1: parse the refseq mappings from reads
		parseReadNamesForRefSeq(mateClusterReadNames, readnameMap);

		//Step 2: parse the repmask mappings from reads
		parseReadNamesForRepMask(mateClusterReadNames, readnameMap, true);
		parseReadNamesForRepMask(splitClusterReadNames, readnameMap, false);

		//Step 3: parse blacklist mappings
		parseReadNamesForBlackList(mateClusterReadNames, readnameMap);
		parseReadNamesForBlackList(splitClusterReadNames, readnameMap);

		//Step 4: parse selfchain mappings
		parseReadNamesForSelfChain(mateClusterReadNames, readnameMap);

		//Step 5: parse inferred insert sizes
		parseReadNamesForInsertSize(mateClusterReadNames, readnameMap);
	}

	protected void parseReadNamesForInsertSize(List<String> readNames,
			Map<String, ReadName> map) {

		for (String readName : readNames){
			ReadName read = map.get(readName);

			if (read != null){
				this.insertSizes.add(read.insertSize);
			}

		}

	}

	protected void parseReadNamesForSelfChain(List<String> readNames,
			Map<String, ReadName> map) {

		for (String readName : readNames){
			ReadName read = map.get(readName);

			if (read != null){
				if (read.selfChainOverlapsMate()){
					this.selfChainOverlap++;
					this.chainScores.add(read.returnScoreOfOverlappingSelfChain());
					if (this.selfChainOverlapString.equals("")){
						this.selfChainOverlapString += read.overlappingSelfChain();
					}else{
						this.selfChainOverlapString += ", " + read.overlappingSelfChain();
					}
				}
			}

		}

	}

	protected void parseReadNamesForBlackList(List<String> readNames,
			Map<String, ReadName> map) {

		for (String readName : readNames){
			ReadName read = map.get(readName);

			if (read != null){
				for (BlacklistAnnotation ann : read.getBlacklist()){
					CollectionUtil.addKeyToCountMap(ann.component, this.blacklistAnchorCounts);
				}
				for (BlacklistAnnotation ann : read.getMateBlacklist()){
					CollectionUtil.addKeyToCountMap(ann.component, this.blacklistMateCounts);
				}
			}
		}

	}

	protected void parseReadNamesForRepMask(List<String> readNames, Map<String, ReadName> map, boolean parseMate) {
		for (String readName : readNames){
			ReadName read = map.get(readName);
			if (read != null){
				for (RepMaskAnnotation repmask : read.getRepMaskAnnotation()){
					String repClass = repmask.repeatClass;
					String repFamily = repmask.repeatFamily;
					CollectionUtil.addKeyToCountMap(repClass, this.repMaskAnchor_counts);
					CollectionUtil.addKeyToCountMap(repFamily, this.repFamilyAnchorCounts);
				}

				if (parseMate){
					for (RepMaskAnnotation repmask : read.getMateRepMaskAnnotation()){
						String repClass = repmask.repeatClass;
						String repFamily = repmask.repeatFamily;
						CollectionUtil.addKeyToCountMap(repClass, this.repMaskMate_counts);
						CollectionUtil.addKeyToCountMap(repFamily, this.repFamilyMateCounts);
					}
				}

			}
		}

	}

	protected void parseReadNamesForRefSeq(List<String> readNames, Map<String, ReadName> map) {
		for (String readName : readNames){
			ReadName read = map.get(readName);
			Set<String> tempSymbols = new HashSet<String>();
			if (read != null){
				for (RefGeneAnnotation gene : read.getMateRefGeneAnnotation()){
					String symbol = gene.geneSymbol;
					if (!tempSymbols.add(symbol)){
						continue; // gene has already be counted, skip to the next one
					}
					CollectionUtil.addKeyToCountMap(symbol, this.refseqMateCounts); //keep track of the gene symbol counts

					if (readName.startsWith(SAMDefinitions.UNIQUE_UNIQUE_MAPPING)){
						CollectionUtil.addKeyToCountMap(symbol, this.refSeqMateUUCounts ); // keep track of gene symbols for uu reads
					}

				}
				for (RefGeneAnnotation gene : read.getRefGeneAnnotation()){
					String symbol = gene.geneSymbol;
					this.refSeqAnchorLocations.add(symbol);
				}
			}
		}
	}

	public List<String> getReadNamesFromMateClusters(){

		List<String> readNames = new ArrayList<String>();

		if(this.hasLeftMateCluster()){
			String names = this.left_cluster_sam_record.getAttribute(MobileDefinitions.SAM_TAG_READNAMES).toString();
			List<String> leftReads = Arrays.asList(names.split(",", -1));
			readNames.addAll(leftReads);
		}
		if(this.hasRightMateCluster()){
			String names = this.right_cluster_sam_record.getAttribute(MobileDefinitions.SAM_TAG_READNAMES).toString();
			List<String> rightReads = Arrays.asList(names.split(",", -1));
			readNames.addAll(rightReads);
		}

		return readNames;


	}

	public List<String> getReadNamesFromSplitClusters(){
		List<String> readNames = new ArrayList<String>();
		if (this.samrecord_splitcluster != null){
			String names = this.samrecord_splitcluster.getAttribute(MobileDefinitions.SAM_TAG_READNAMES).toString();
			readNames = Arrays.asList(names.split(",", -1));
		}
		return readNames;
	}

	public void parseSAMRecordCluster(SAMRecord cluster) throws IllegalSAMPairException {
		Boolean isSplitCluster;
		String mobileHit;
		mobileHit = cluster.getAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT).toString();
		this.parseSampleCountTag(cluster.getAttribute(MobileDefinitions.SAM_TAG_SAMPLECOUNT).toString());


		for (String mobileMapping : mobileHit.split(";", -1)){
			if (!mobileHit.equals("")){
				this.mobile_mappings.add(mobileMapping);
			}
		}

		if (this.original_reference.equals("")){
			this.original_reference = cluster.getReferenceName();
		}else if(!this.original_reference.equals(cluster.getReferenceName())){
			throw new IllegalSAMPairException("Could not add cluster to MobilePrediction. " +
					"Clusters are not originating from same reference: " + cluster.getReferenceName()
					+ "vs" + this.original_reference );
		}



		isSplitCluster = Boolean.parseBoolean(cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_CLUSTER).toString());
		if (isSplitCluster){
			parseSplitCluster(cluster);
		}else{
			parseMateCluster(cluster);
		}
	}

	protected void parseSampleCountTag(String tag){

		String[] pairs = tag.split(",", -1);
		String sample;
		int count;
		int oldCount;


		if (!"".equals(pairs[0])){
			for (String pair : pairs){
				String[] sampleCount = pair.split("=");
				sample = sampleCount[0].trim();
				count = Integer.parseInt(sampleCount[1]);
				this.sample_names.add(sample);

				if (this.sample_counts.containsKey(sample)){
					oldCount = this.sample_counts.get(sample);
					this.sample_counts.put(sample, oldCount + count);
				}else{
					this.sample_counts.put(sample, count);
				}
			}
		}
	}

	protected void parseMateCluster(SAMRecord cluster){
		String clusterHits;
		String clusterLength;

		this.discordant_read_names.add(cluster.getReadName());

		clusterHits = cluster.getAttribute(MobileDefinitions.SAM_TAG_CLUSTER_HITS).toString();
		clusterLength = cluster.getAttribute(MobileDefinitions.SAM_TAG_CLUSTER_LENGTH).toString();

		this.unique_hits += Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_UNIQUE_HITS).toString());
		this.multiple_hits += Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_MULTIPLE_HITS).toString());
		this.unmapped_hits += Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_UNMAPPED_HITS).toString());

		this.mate_mobile_mapping = cluster.getAttribute(MobileDefinitions.SAM_TAG_MOBILE_HIT).toString();

		if(!cluster.getReadNegativeStrandFlag()){
			this.left_cluster_sam_record = cluster;
			this.left_cluster_length = Integer.parseInt(clusterLength);
			this.left_mate_hits = Integer.parseInt(clusterHits);
			this.left_mate_cluster_border = cluster.getAlignmentEnd();
			this.left_mate_median_mapq = cluster.getMappingQuality();
		}else{
			this.right_cluster_sam_record = cluster;
			this.right_cluster_length = Integer.parseInt(clusterLength);
			this.right_mate_hits = Integer.parseInt(clusterHits);
			this.right_mate_cluster_border = cluster.getAlignmentStart() - 1;
			this.right_mate_median_mapq = cluster.getMappingQuality();
		}
	}


	protected void parseSplitCluster(SAMRecord cluster){
		String leftClippedHits;
		String rightClippedHits;
		String leftClippedMedianEnd;
		String rightClippedMedianEnd;

		this.split_read_names.add(cluster.getReadName());

		leftClippedHits = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_HITS).toString();
		this.right_aligned_split_hits = Integer.parseInt(leftClippedHits); //left clipped reads are aligned to right of MEI
		leftClippedMedianEnd = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_MEDIAN_END).toString();
		if (Integer.parseInt(leftClippedMedianEnd) != 0){
			this.right_aligned_split_border = Integer.parseInt(leftClippedMedianEnd) - 1;
		}

		rightClippedHits = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_HITS).toString();
		this.left_aligned_split_hits = Integer.parseInt(rightClippedHits); //right clipped reads are aligned to left of MEI
		rightClippedMedianEnd = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_MEDIAN_END).toString();
		this.left_aligned_split_border = Integer.parseInt(rightClippedMedianEnd);

		this.leftclipped_max_distance = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_MAX_DISTANCE).toString();
		this.rightclipped_max_distance = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_MAX_DISTANCE).toString();
		this.leftclipped_same_fraction = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_FRAC_SAME_DISTANCE).toString();
		this.rightclipped_same_fraction = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_FRAC_SAME_DISTANCE).toString();
		this.clipped_avg_qual = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_CLIPPED_AVG_QUAL).toString();
		this.clipped_avg_len = cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_AVG_CLIPPED_LEN).toString();

		//right and left are switches before and after assignment (=) operator, because rightclipped reads are left aligned!
		this.left_aligned_polyA_hits = Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_POLYA).toString());
		this.left_aligned_polyT_hits = Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_POLYT).toString());
		this.right_aligned_polyA_hits = Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_POLYA).toString());
		this.right_aligned_polyT_hits = Integer.parseInt(cluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_POLYT).toString());
		this.samrecord_splitcluster = cluster;
	}

	protected void createFeatures(){
		for (String feature : HEADER){
			if(feature.equals(COLUMN_BORDER3)){
				features.put(feature, Integer.toString(this.getRightPredictionBorder()));
			}else if (feature.equals(COLUMN_BORDER5)) {
				features.put(feature, Integer.toString(this.getLeftPredictionBorder()));
			} else if (feature.equals(COLUMN_ENDBORDER3)){
				features.put(feature, Integer.toString(this.getEndRightPredictionBorder()));
			}else if (feature.equals(COLUMN_ENDBORDER5)){
				features.put(feature, Integer.toString(this.getEndLeftPredictionBorder()));
			}else if (feature.equals(COLUMN_CLUSTER3_HITS)){
				String cluster3Hits = (this.getRightTotalHits()== 0) ? "NA" : Integer.toString(this.getRightTotalHits());
				features.put(feature, cluster3Hits);
			}else if (feature.equals(COLUMN_MERGED)){
				features.put(feature, Boolean.toString(this.merged));
			}
			else if (feature.equals(COLUMN_CLUSTER3_LEN)){
				String cluster3Len = (this.right_cluster_length == 0) ? "NA" : Integer.toString(this.right_cluster_length);
				features.put(feature, cluster3Len);
			}else if (feature.equals(COLUMN_CLUSTER5_HITS)){
				String cluster5Hits = (this.getLeftTotalHits() == 0) ? "NA" : Integer.toString(this.getLeftTotalHits());
				features.put(feature, cluster5Hits);
			}else if (feature.equals(COLUMN_CLUSTER5_LEN)){
				String cluster5Len = (this.left_cluster_length == 0) ? "NA" : Integer.toString(this.left_cluster_length);
				features.put(feature, cluster5Len);
			}else if (feature.equals(COLUMN_INSERTPOINT)){
				features.put(COLUMN_INSERTPOINT, Integer.toString(this.getInsertionEstimate()));
			}else if (feature.equals(COLUMN_ENDINSERTPOINT)){
				features.put(COLUMN_ENDINSERTPOINT, Integer.toString(this.getEndInsertionEstimate()));
			}else if (feature.equals(COLUMN_MOBILE)){
				features.put(feature, CollectionUtil.toString(this.mobile_mappings));
			}else if (feature.equals(COLUMN_REFERENCE)){
				features.put(feature, this.original_reference);
			}else if (feature.equals(COLUMN_SPLIT5_HITS)){
				features.put(feature, Integer.toString(this.left_aligned_split_hits));
			}else if (feature.equals(COLUMN_SPLIT3_HITS)){
				features.put(feature, Integer.toString(this.right_aligned_split_hits));
			}else if (feature.equals(COLUMN_UNIQUE_HITS)){
				features.put(feature, Integer.toString(this.unique_hits));
			}else if (feature.equals(COLUMN_MULTIPLE_HITS)){
				features.put(feature, Integer.toString(this.multiple_hits));
			}else if (feature.equals(COLUMN_UNMAPPED_HITS)){
				features.put(feature, Integer.toString(this.unmapped_hits));
			}else if (feature.equals(COLUMN_LEFTCLIPPED_MAX_DISTANCE)){
				features.put(feature, this.leftclipped_max_distance);
			}else if (feature.equals(COLUMN_RIGHTCLIPPED_MAX_DISTANCE)){
				features.put(feature,  this.rightclipped_max_distance);
			}else if (feature.equals(COLUMN_LEFTCLIPPED_FRAC_DISTANCE)){
				features.put(feature, this.leftclipped_same_fraction);
			}else if (feature.equals(COLUMN_RIGHTCLIPPED_FRAC_DISTANCE)){
				features.put(feature, this.rightclipped_same_fraction);
			}else if (feature.equals(COLUMN_CLIPPED_AVG_QUAL)){
				features.put(feature, this.clipped_avg_qual);
			}else if (feature.equals(COLUMN_TSD)){
				features.put(feature, this.hasTSD());
			}else if (feature.equals(COLUMN_TSD_LEN)){
				features.put(feature, Integer.toString(this.getTSDlength()));
			}else if (feature.equals(COLUMN_TSD_SEQ)){
				NucleotideSequence TSDseq = getTSDseq();
				if(TSDseq != null)
					features.put(COLUMN_TSD_SEQ, TSDseq.getSequence());
				else features.put(feature, ".");
			}else if (feature.equals(COLUMN_AVG_CLIPPED_LEN)){
				features.put(feature, this.clipped_avg_len);
			}else if (feature.equals(COLUMN_SAMPLE)){
				String sampleNames = this.sample_names.toString();
				sampleNames = sampleNames.substring(1, sampleNames.length() - 1);
				features.put(feature, sampleNames);
			}else if (feature.equals(COLUMN_POLYA3_HITS)){
				features.put(feature, Integer.toString(this.right_aligned_polyA_hits));
			}else if (feature.equals(COLUMN_POLYT3_HITS)){
				features.put(feature, Integer.toString(this.right_aligned_polyT_hits));
			}else if (feature.equals(COLUMN_POLYA5_HITS)){
				features.put(feature, Integer.toString(this.left_aligned_polyA_hits));
			}else if (feature.equals(COLUMN_POLYT5_HITS)){
				features.put(feature, Integer.toString(this.left_aligned_polyT_hits));
			}else if (feature.equals(COLUMN_SAMPLE_COUNT)){
				String sampleCount = this.sample_counts.toString();
				sampleCount = sampleCount.substring(1, sampleCount.length() - 1);
				features.put(feature, sampleCount);
			} else if (feature.equals(COLUMN_NON_SUPPORTING_SAMPLE_COUNT)){
				String sampleNonSupportingCount = this.non_supporting_sample_counts.toString();
				sampleNonSupportingCount = sampleNonSupportingCount.substring(1, sampleNonSupportingCount.length() - 1);
				features.put(feature, sampleNonSupportingCount);
			} else if(feature.equals(COLUMN_VAF)) {
				String sampleVAFstring = this.sample_vafs.toString();
				features.put(feature, sampleVAFstring.substring(1, sampleVAFstring.length() - 1));
			} else if(feature.equals(COLUMN_SOMATIC)){
				if(this.output_somatic)
					features.put(feature, this.somatic);
			}
		}
	}

	public String getOriginalReference(){
		return this.original_reference;
	}

	public Map<String, Integer> getSampleCounts(){
		return this.sample_counts;
	}

	public Set<String> getSampleNames(){
		return this.sample_names;
	}

	public String getSamplesNamesAsString(){
		String sampleNames = this.sample_names.toString();
		sampleNames = sampleNames.substring(1, sampleNames.length() - 1);
		return sampleNames;
	}

	public String getSampleCountsAsString(){
		String sampleCount = this.sample_counts.toString();
		sampleCount = sampleCount.substring(1, sampleCount.length() - 1);
		return sampleCount;
	}

	public int getLeftMateClusterLength(){
		return this.left_cluster_length;
	}

	public int getRightMateClusterLength(){
		return this.right_cluster_length;
	}

	public int getRightMateClusterBorder(){ return this.right_mate_cluster_border;}

	public int getLeftMateClusterBorder(){
		return this.left_mate_cluster_border;
	}

	public String hasTSD(){
		if (this.hasLeftAlignedSplitCluster() && this.hasRightAlignedSplitCluster()){
			if (this.right_aligned_split_border < this.left_aligned_split_border){
				return DUPLICATION;
			}else if(this.right_aligned_split_border > this.left_aligned_split_border){
				return DELETION;
			}else{
				return NO_TSD;
			}

		}else if (this.hasLeftAlignedSplitCluster() && this.hasRightMateCluster()){
			if (this.right_mate_cluster_border < this.left_aligned_split_border){
				return DUPLICATION;
			}
		}else if (this.hasRightAlignedSplitCluster() && this.hasLeftMateCluster()){
			if (this.right_aligned_split_border < this.left_mate_cluster_border){
				return DUPLICATION;
			}
		}else if (this.hasLeftMateCluster() && this.hasRightMateCluster()){
			if (this.right_mate_cluster_border < this.left_mate_cluster_border){
				return DUPLICATION;
			}
		}

		return UNKNOWN;
	}

	protected int getHighestBorder() {
		Vector<Integer> borders = new Vector<Integer>();

		if (this.hasLeftAlignedSplitCluster() && this.hasRightAlignedSplitCluster()){
			borders.add(this.left_aligned_split_border);
			borders.add(this.right_aligned_split_border);
			Collections.sort(borders);
			return borders.lastElement();
		}
		//else
		borders.add(this.left_aligned_split_border);
		borders.add(this.right_aligned_split_border);
		borders.add(this.left_mate_cluster_border);
		borders.add(this.right_mate_cluster_border);
		Collections.sort(borders);
		return borders.lastElement();
	}

	public String getChromosome(){
		return this.original_reference;
	}

	protected int getLowestBorder(){
		Vector<Integer> borders = new Vector<Integer>();
		
		if (this.hasLeftAlignedSplitCluster() && this.hasRightAlignedSplitCluster()){
			borders.add(this.left_aligned_split_border);
			borders.add(this.right_aligned_split_border);
			Collections.sort(borders);
			return borders.firstElement();
		}
		if (this.left_aligned_split_border != 0){
			borders.add(this.left_aligned_split_border);
		}
		if (this.right_aligned_split_border != 0){
			borders.add(this.right_aligned_split_border);
		}
		if (this.left_mate_cluster_border != 0){
			borders.add(this.left_mate_cluster_border);
		}
		if (this.right_mate_cluster_border != 0){
			borders.add(this.right_mate_cluster_border);
		}
		Collections.sort(borders);
		return borders.firstElement();
	}

	public int getLeftPredictionBorder(int extraWindow){
		return getLeftPredictionBorder() - extraWindow;
	}

	public int getRightPredictionBorder(int extraWindow){
		return getRightPredictionBorder() + extraWindow;
	}

	public int getEndLeftPredictionBorder(int extraWindow){
		return getEndLeftPredictionBorder() - extraWindow;
	}
	public int getEndRightPredictionBorder(int extraWindow){
		return getEndRightPredictionBorder() + extraWindow;
	}

	public Set<String> getDiscordantClusterIds(){
		return this.discordant_read_names;
	}
	
	public Set<String> getSplitClusterIds(){
		return this.split_read_names;
	}
	
	public int getLeftTotalHits(){
		return (this.left_aligned_split_hits + this.left_mate_hits);
	}
	
	public int getRightTotalHits(){
		return (this.right_aligned_split_hits + this.right_mate_hits);
	}
	
	public int getLeftPolyAhits(){
		return this.left_aligned_polyA_hits;
	}
	
	public int getLeftPolyThits(){
		return this.left_aligned_polyT_hits;
	}
	
	public int getRightPolyAhits(){
		return this.right_aligned_polyA_hits;
	}
	
	public int getRightPolyThits(){
		return this.right_aligned_polyT_hits;
	}
	
	public boolean hasOneOrZeroHomoPolymerMappings(){
		if (this.left_aligned_polyA_hits + this.left_aligned_polyT_hits == 0){
			return true;
		}else if (this.right_aligned_polyA_hits + this.right_aligned_polyT_hits == 0){
			return true;
		}
		return false;
	}
	
	
	public boolean hasLeftMateCluster(){
		return (this.left_mate_hits > 0);
	}
	
	public boolean hasLeftAlignedSplitCluster(){
		return (this.left_aligned_split_hits > 0);
	}
	
	
	public boolean hasRightAlignedSplitCluster(){
		return (this.right_aligned_split_hits > 0);
	}
	
	public boolean hasRightMateCluster(){
		return (this.right_mate_hits > 0);
	}
	
	public SimpleRegion predictionWindowToRegion(){
		
		int leftBorder = this.getLeftPredictionBorder();
		int rightBorder = this.getEndRightPredictionBorder();
		
		//TODO: this correction should be moved into .getLeftPrediction and .getRightPrediction
		if (leftBorder < 1){
			leftBorder = 1;
		}
		
		if (rightBorder < 2){
			rightBorder = 2;
		}
		
		return new SimpleRegion(this.getOriginalReference(), leftBorder, rightBorder);
	}
	
	public SAMRecord getSplitSAMRecord(){
		return this.samrecord_splitcluster;
	}
	
	public String getMateMobileMapping(){
		return this.mate_mobile_mapping;
	}
	
	public Set<String> getMobileMappings(){
		return this.mobile_mappings;
	}

	public static String getHeader(){
		return getHeader(false);
	}

	public static String getHeader(boolean outputSomatic){
	StringBuilder sb = new StringBuilder();
	for (int i = 0; i < HEADER.length; i++){
		if(!HEADER[i].equals(COLUMN_SOMATIC) || (outputSomatic && HEADER[i].equals(COLUMN_SOMATIC))){
			sb.append(HEADER[i]);
			sb.append("\t");
		}
	}
	return sb.toString().trim();
}
	
	public String toString(){
		createFeatures();
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < HEADER.length; i++){
			if(!HEADER[i].equals(COLUMN_SOMATIC) || (output_somatic && HEADER[i].equals(COLUMN_SOMATIC))){
				sb.append(features.get(HEADER[i]));
				sb.append("\t");
			}
		}
		return sb.toString().trim();
	}
	
	public List<Integer> getInsertsFromSupportingReads(boolean doNotAddZeros){
	List<Integer> inserts = new ArrayList<Integer>();
		
		//remove 0s from insertsizes
		for (int i : this.insertSizes){
			if (i != 0 && doNotAddZeros){
				inserts.add(i);
			}
		}
		return inserts;
	}
	
	public List<Integer> getInsertsFromSupportingReads(){
		return this.getInsertsFromSupportingReads(true);
	}

	public void setNonSupportingCount(String sampleName, int count){ this.non_supporting_sample_counts.put(sampleName, count); }

	public Map<String, Integer> getNonSupportingCounts(String sampleName, double count){ return this.non_supporting_sample_counts; }

	public void setVAF(String sampleName, double VAF){
		this.sample_vafs.put(sampleName, VAF);
	}

	public Map<String, Double> getVAFs(String sampleName){
		return this.sample_vafs;
	}
	
	public void setReferenceGenome(ReferenceGenome referenceGenome){
		this.reference_genome = referenceGenome;
	}

	public String toGripsString(){
		Set<String> overlappingRepClass = CollectionUtil.returnOverlappingKeysInMap(this.repMaskAnchor_counts, this.repMaskMate_counts);
		Set<String> overlappingRepFamily = CollectionUtil.returnOverlappingKeysInMap(this.repFamilyAnchorCounts, this.repFamilyMateCounts);
		
		List<Integer> insertsNoZeros = this.getInsertsFromSupportingReads();
		
		return this.getChromosome() + "\t" + this.getLeftPredictionBorder() + "\t" + this.getRightPredictionBorder() +
				"\t" + Integer.toString(this.getLeftTotalHits() + this.getRightTotalHits()) + "\t" + CollectionUtil.toString(this.mobile_mappings) + "\t" + this.getInsertionEstimate() + "\t"
				+ this.getSamplesNamesAsString() + "\t" + this.getSampleCountsAsString() + "\t" + this.left_cluster_length + "\t" + this.right_cluster_length + "\t" 
				+ Integer.toString(this.left_aligned_split_hits + this.right_aligned_split_hits) + "\t" + this.getLeftTotalHits() + "\t" + this.getRightTotalHits() + "\t" + CollectionUtil.mapToString(this.refseqMateCounts) + "\t" + CollectionUtil.mapToString(this.repMaskAnchor_counts) +
				"\t" +  CollectionUtil.mapToString(this.repMaskMate_counts) + "\t" + this.sameRefSeqMappingAsPrediction() + "\t" +
				overlappingRepClass.toString().substring(1, overlappingRepClass.toString().length() - 1) + "\t"
				+ sameAnchorLocationAsPrediction() + "\t" + CollectionUtil.mapToString(this.blacklistAnchorCounts) + "\t" +
				CollectionUtil.mapToString(this.blacklistMateCounts) + "\t" + CollectionUtil.mapToString(this.repFamilyAnchorCounts) + "\t" + 
				CollectionUtil.mapToString(this.repFamilyMateCounts) + "\t" + overlappingRepFamily.toString().substring(1, overlappingRepFamily.toString().length() - 1) +
				"\t" + this.unique_hits + "\t" + this.multiple_hits + "\t" + this.unmapped_hits + "\t" + this.selfChainOverlap + "\t" +
				this.selfChainOverlapString + "\t" + this.chainScores.toString() + "\t" + MathFunction.getMedianFromDoubles(this.chainScores) +
				"\t" + this.insertSizes.toString() + "\t" + MathFunction.getMedianFromIntegers(insertsNoZeros) + "\t" +
				Boolean.toString(this.sameRefSeqMappingAsPrediction() || this.hasOnlyDiscordantUXReads(MIN_DISCORDANT_UX)) + "\t" + Integer.toString(this.left_aligned_polyA_hits + this.left_aligned_polyT_hits) +
				"\t" + Integer.toString(this.right_aligned_polyA_hits + this.right_aligned_polyT_hits) + "\t" + this.hasTSD() + "\t" + this.left_mate_median_mapq
				+ "\t" + this.right_mate_median_mapq;
	}
	
	public boolean hasOnlyDiscordantUXReads(int uxThreshold){
		return (this.multiple_hits == 0 && this.unique_hits == 0 && this.unmapped_hits >= uxThreshold);
	}
	
	public boolean hasOnlyPercentageDiscordantUXReads(double uxThreshold){
		int discordantHits = this.multiple_hits + this.unique_hits + this.unmapped_hits;
		double perc = (double) this.unmapped_hits / (double) discordantHits * 100;
		return perc >= uxThreshold;
	}
	
	
	public boolean gripNeedsFiltering(){
		
		List<Integer> inserts = this.getInsertsFromSupportingReads();
		
		//RefSeq prediction is same as where anchor is mapped
		if (this.sameAnchorLocationAsPrediction()){
			return true;
		}
		//If there is overlap with blacklist
		if (! this.blacklistAnchorCounts.isEmpty() || ! this.blacklistMateCounts.isEmpty()){
			return true;
		}
		//if there is selfchaining
		if (! this.chainScores.isEmpty()){
			return true;
		}
		// if there is overlapping rep family in anchor and mate
		if (! CollectionUtil.returnOverlappingKeysInMap(this.repFamilyAnchorCounts, this.repFamilyMateCounts).isEmpty()){
			return true;
		}
		
		// if it overlaps centromeres
		if (this.repFamilyAnchorCounts.containsKey("centr") || this.repFamilyMateCounts.containsKey("centr")){
			return true;
		}
		
		// if it is not mapping to refseq in original bam file and if this can not be explained by present UX reads
		if (! this.sameRefSeqMappingAsPrediction() && ! this.hasOnlyDiscordantUXReads(MIN_DISCORDANT_UX)){
			return true;
		}
		
//		if (! this.sameRefSeqMappingAsPrediction() && ! this.hasOnlyPercentageDiscordantUXReads(62)){
//			return true;
//		}
		
		if (MathFunction.getMedianFromIntegers(inserts) != null && MathFunction.getMedianFromIntegers(inserts) < MIN_INSERT_SIZE){
			return true;
		}
		
		
		
		return false;
	}
	
	public String toGripsHeader(){
		return "Chr\tStart\tEnd\tReads\tPrediction\tInsertion\tSample\tSample counts\tLeftClusterLength\tRightClusterLength\tTotalSplitHits\tLeftTotalHits\tRightTotalHits\tRefSeqMateCount\tRepMaskAnchorCount\tRepMaskMateCount\tSameRefSeq\t" +
				"OverlappingRepeatClass\tsameAnchorLocationAsPred\tBlackListAnchorCount\tBlackListMateCount\tRepFamilyAnchorCount\t" +
				"RepFamilyMateCount\tOverlappingRepFamily\tUU\tUM\tUX\tSelfChain\tSelfChainDetailed\tSelfChainScores\tSelfChainScoresMedian\t" +
				"insert sizes\tmedian insert\tSameRefSeqOrUX\tLeftAlignedPoly\tRightAlignedPoly\tTSD\tleftMateMedianMAPQ\trightMateMedianMAPQ";
	}
	
	public boolean sameRefSeqUUMappingAsPrediction(){
		for (String prediction : this.mobile_mappings){
			if (this.refSeqMateUUCounts.keySet().contains(prediction)){
				return true;
			}
		}
		
		return false;
	}
	
	public boolean sameRefSeqMappingAsPrediction(){
		
		for (String prediction : this.mobile_mappings){
			if(this.refseqMateCounts.keySet().contains(prediction)){
				return true;
			}
		}
		
		return false;

	}
	
	public boolean sameAnchorLocationAsPrediction(){
		for (String prediction : this.mobile_mappings){
			if(this.refSeqAnchorLocations.contains(prediction)){
				return true;
			}
		}
		
		return false;

	}
	

	

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime
				* result
				+ ((discordant_read_names == null) ? 0 : discordant_read_names
						.hashCode());
		result = prime * result + left_aligned_split_border;
		result = prime * result + left_aligned_split_hits;
		result = prime * result + left_cluster_length;
		result = prime * result + left_mate_cluster_border;
		result = prime * result + left_mate_hits;
		result = prime * result
				+ ((mobile_mappings == null) ? 0 : mobile_mappings.hashCode());
		result = prime
				* result
				+ ((original_reference == null) ? 0 : original_reference
						.hashCode());
		result = prime * result + right_aligned_split_border;
		result = prime * result + right_aligned_split_hits;
		result = prime * result + right_cluster_length;
		result = prime * result + right_mate_cluster_border;
		result = prime * result + right_mate_hits;
		result = prime * result
				+ ((sample_names == null) ? 0 : sample_names.hashCode());
		result = prime
				* result
				+ ((split_read_names == null) ? 0 : split_read_names.hashCode());
		return result;
	}
	
	public boolean onlyLeftDiscordantReads(){
		return (this.left_mate_hits > 0 && this.right_mate_hits == 0);
	}
	
	public boolean onlyRightDiscordantReads(){
		return (this.right_mate_hits > 0 && this.left_mate_hits == 0);
	}


	public int getLeftClusterBorder(){
		if(hasLeftAlignedSplitCluster())
			return left_aligned_split_border;
		else if (hasLeftMateCluster())
			return left_mate_cluster_border;
		else return 0;
	}
	public int getRightClusterBorder(){
		if(hasRightAlignedSplitCluster())
			return right_aligned_split_border;
		else if (hasRightMateCluster())
			return right_mate_cluster_border;
		else return 0;
	}

	public boolean hasLeftAndRightCluster(){
		return (hasLeftAlignedSplitCluster() || hasLeftMateCluster())
				&& (hasRightAlignedSplitCluster() || hasRightMateCluster());
	}

	public int getTSDlength(){
		if(!Objects.equals(hasTSD(), UNKNOWN))
			return Math.abs(getLeftClusterBorder() - getRightClusterBorder());
		else
			return -1;
	}

	public NucleotideSequence getTSDseq() {
		if(!hasTSD().equals(UNKNOWN) && !hasTSD().equals(NO_TSD) && this.reference_genome != null)
			try {
				return reference_genome.getSubSequenceAt(getChromosome(), getInsertionEstimate() + 1, getEndInsertionEstimate());
			} catch (InvalidNucleotideSequenceException e) {
				return null;
			}
		else
			return null;
	}

	private int getSingleClusterInsertionEstimate(){

		//When there are both left and right clusters, raise an exception
		if(hasLeftAndRightCluster())
			throw new IllegalStateException("There are no single clusters on only the left and right side.");

		//Just return the border for split clusters
		if(this.hasLeftAlignedSplitCluster())
			return left_aligned_split_border;
		else if(this.hasRightAlignedSplitCluster())
			return right_aligned_split_border;

		//Correct the border for mate clusters if necessary
		else if(hasLeftMateCluster())
			if(this.left_cluster_length >= median_fragment_length + sd_fragment_length)
				return left_mate_cluster_border;
			else return (left_mate_cluster_border +
					left_mate_cluster_border + (median_fragment_length - left_cluster_length) + sd_fragment_length) / 2;
		else if(hasRightMateCluster())
			if(right_cluster_length >= median_fragment_length + sd_fragment_length)
				return right_mate_cluster_border;
			else return (right_mate_cluster_border + (right_mate_cluster_border -
					(median_fragment_length - right_cluster_length) - sd_fragment_length)) / 2;
		else return 0;
	}


//		insertionEstimate = (this.left_mate_cluster_border +
//		this.left_mate_cluster_border + (median_fragment_length - this.left_cluster_length) + sd_fragment_length) / 2;
//		}
//		}else if(hasRightMateCluster() && !hasLeftMateCluster()){
//		if (this.right_cluster_length >= median_fragment_length + sd_fragment_length){
//		insertionEstimate = this.right_mate_cluster_border;
//		}else{
//		insertionEstimate = (this.right_mate_cluster_border + (this.right_mate_cluster_border -
//		(median_fragment_length - this.right_cluster_length) - sd_fragment_length)) / 2;
//		}

	public int getInsertionEstimate(){
		int insertionEstimate = 0;

		//Left+right clusters: When a target site duplication is present, return the right border
		//Otherwise return the left split border, right split border or for two mate clusters their middle
		if(hasLeftAndRightCluster())
			if(hasTSD().equals(DUPLICATION))
				insertionEstimate = getRightClusterBorder();
			else if(hasLeftAlignedSplitCluster())
				insertionEstimate = left_aligned_split_border;
			else if(hasRightAlignedSplitCluster())
				insertionEstimate = right_aligned_split_border;
			else
				insertionEstimate = Math.round((this.left_mate_cluster_border + this.right_mate_cluster_border) / 2);

		//Single clusters: Return only their border
		else
			insertionEstimate = getSingleClusterInsertionEstimate();

		return insertionEstimate;
	}

	public int getEndInsertionEstimate(){
		int endInsertionEstimate = 0;

		//Left+right clusters: When a target site duplication is present, return the left border
		//Otherwise return the right split border, left split border or for two mate clusters their middle
		if(hasLeftAndRightCluster())
			if(hasTSD().equals(DUPLICATION))
				endInsertionEstimate = getLeftClusterBorder();
			else if(hasRightAlignedSplitCluster())
				endInsertionEstimate = right_aligned_split_border;
			else if(hasLeftAlignedSplitCluster())
				endInsertionEstimate = left_aligned_split_border;
			else
				endInsertionEstimate = Math.round((this.left_mate_cluster_border + this.right_mate_cluster_border) / 2);

			//Single clusters: Return only their border
		else
			endInsertionEstimate = getSingleClusterInsertionEstimate();

		return endInsertionEstimate;
	}

	public int getLeftPredictionBorder(){
		if(max_expected_cluster_size == 1333)
			logger.error("ERROR!");
		int leftPredictionBorder = 0;

		//For a left+right split cluster, return the insertion point
		if(hasLeftAlignedSplitCluster() && hasRightAlignedSplitCluster())
			leftPredictionBorder = getInsertionEstimate();

		//For a left split + right mate cluster with duplication,
		//the border is the mate cluster border increased up to the expected target site duplication length
		//For a left mate + right split with duplication,
		//the border is the that of the right split cluster
		//In case there is no duplication, the lowest coordinate is taken
		else if(hasLeftAlignedSplitCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				int TSDlength = getTSDlength();
				if(TSDlength < EXPECTED_DUPLICATION_LENGTH)
					leftPredictionBorder = right_mate_cluster_border - (EXPECTED_DUPLICATION_LENGTH - TSDlength);
				else
					leftPredictionBorder = right_mate_cluster_border;
			} else
				leftPredictionBorder = getLowestBorder();
		} else if(hasLeftMateCluster() && hasRightAlignedSplitCluster()){
			if(hasTSD().equals(DUPLICATION)){
				leftPredictionBorder = right_aligned_split_border;
			} else
				leftPredictionBorder = getLowestBorder();
		}

		//For two mate clusters with duplication,
		//the border is increased up to half of expected target site duplication length
		//Otherwise the lowest coordinate is taken
		else if(hasLeftMateCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				int TSDlength = getTSDlength();
				if(TSDlength < EXPECTED_DUPLICATION_LENGTH)
					leftPredictionBorder = right_mate_cluster_border - (EXPECTED_DUPLICATION_LENGTH - TSDlength) / 2;
				else
					leftPredictionBorder = right_mate_cluster_border;
			} else
				leftPredictionBorder = getLowestBorder();
		}

		//For single clusters substract the fuzziness
		//unless the size of the right mate cluster is unexpectedly small. Then increase it to the expected size
		else if(hasRightAlignedSplitCluster()){
			leftPredictionBorder =  this.right_aligned_split_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasLeftAlignedSplitCluster()){
			leftPredictionBorder =  this.left_aligned_split_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasLeftMateCluster()){
			leftPredictionBorder =  this.left_mate_cluster_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else{
			if (this.right_cluster_length >= this.max_expected_cluster_size - SINGLE_CLUSTER_BORDER_FUZZINESS){
				leftPredictionBorder = this.right_mate_cluster_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
			}
			else{
				leftPredictionBorder = (this.right_mate_cluster_border - (this.max_expected_cluster_size - this.right_cluster_length));
			}
		}

		return leftPredictionBorder;
	}

	public int getRightPredictionBorder(){
		int rightPredictionBorder = 0;

		//For a left+right split cluster, return the insertion point
		if(hasLeftAlignedSplitCluster() && hasRightAlignedSplitCluster())
			rightPredictionBorder = getInsertionEstimate();

			//For a split + mate cluster with duplication,
			//the border is the mate cluster border increased up to the expected target site duplication length
			//Otherwise the highest coordinate is taken
		else if(hasLeftAlignedSplitCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				rightPredictionBorder = right_mate_cluster_border;
			} else
				rightPredictionBorder = getHighestBorder();
		} else if(hasLeftMateCluster() && hasRightAlignedSplitCluster()){
			if(hasTSD().equals(DUPLICATION)){
				rightPredictionBorder = right_aligned_split_border;
			} else
				rightPredictionBorder = getHighestBorder();
		}

		//For two mate clusters with duplication,
		//the border is the mate cluster border
		//Otherwise the lowest coordinate is taken
		else if(hasLeftMateCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				rightPredictionBorder = right_mate_cluster_border;
			} else
				rightPredictionBorder = getHighestBorder();
		}

		//For single clusters add the fuzziness
		//unless the size of the right mate cluster is unexpectedly small. Then increase it to the expected size
		else if(hasLeftAlignedSplitCluster()){
			rightPredictionBorder = this.left_aligned_split_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasRightAlignedSplitCluster()){
			rightPredictionBorder = this.right_aligned_split_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasRightMateCluster()){
			rightPredictionBorder = this.right_mate_cluster_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}
		else{
			if(this.left_cluster_length >= this.max_expected_cluster_size - SINGLE_CLUSTER_BORDER_FUZZINESS){
				rightPredictionBorder = (this.left_mate_cluster_border + SINGLE_CLUSTER_BORDER_FUZZINESS);
			}
			else
				rightPredictionBorder = (this.left_mate_cluster_border + this.max_expected_cluster_size - this.left_cluster_length);
		}

		return rightPredictionBorder;
	}

	public int getEndLeftPredictionBorder(){
		int endLeftPredictionBorder = 0;

		//For a left+right split cluster, return the end insertion point
		if(hasLeftAlignedSplitCluster() && hasRightAlignedSplitCluster())
			endLeftPredictionBorder = getEndInsertionEstimate();

		//For a split + mate cluster with duplication,
		//the border is the mate cluster border increased up to the expected target site duplication length
		//Otherwise the lowest coordinate is taken
		else if(hasLeftMateCluster() && hasRightAlignedSplitCluster()){
			if(hasTSD().equals(DUPLICATION)){
				endLeftPredictionBorder = left_mate_cluster_border;
			} else
				endLeftPredictionBorder = getLowestBorder();
		} else if(hasLeftAlignedSplitCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				endLeftPredictionBorder = left_aligned_split_border;
			} else
				endLeftPredictionBorder = getLowestBorder();
		}

		//For two mate clusters with duplication,
		//the border is the mate cluster border
		//Otherwise the highest coordinate is taken
		else if(hasLeftMateCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				endLeftPredictionBorder = left_mate_cluster_border;
			} else
				endLeftPredictionBorder = getLowestBorder();
		}

		//For single clusters substract the fuzziness
		//unless the size of the right mate cluster is unexpectedly small. Then increase it to the expected size
		else if (hasRightAlignedSplitCluster()){
			endLeftPredictionBorder = this.right_aligned_split_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if (hasLeftAlignedSplitCluster()){
			endLeftPredictionBorder = this.left_aligned_split_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if (hasLeftMateCluster()){
			endLeftPredictionBorder = this.left_mate_cluster_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else{
			if (this.right_cluster_length >= this.max_expected_cluster_size - SINGLE_CLUSTER_BORDER_FUZZINESS){
				endLeftPredictionBorder = this.right_mate_cluster_border - SINGLE_CLUSTER_BORDER_FUZZINESS;
			}
			else
				endLeftPredictionBorder = (this.right_mate_cluster_border - (this.max_expected_cluster_size - this.right_cluster_length));
		}

		return endLeftPredictionBorder;
	}

	public int getEndRightPredictionBorder(){
		int endRightPredictionBorder = 0;

		//For a left+right split cluster, return the end insertion point
		if(hasLeftAlignedSplitCluster() && hasRightAlignedSplitCluster())
			endRightPredictionBorder = getEndInsertionEstimate();

			//For a split + mate cluster with duplication,
			//the border is the mate cluster border increased up to the expected target site duplication length
			//Otherwise the highest coordinate is taken
		else if(hasLeftMateCluster() && hasRightAlignedSplitCluster()){
			if(hasTSD().equals(DUPLICATION)){
				int TSDlength = getTSDlength();
				if(TSDlength < EXPECTED_DUPLICATION_LENGTH)
					endRightPredictionBorder = left_mate_cluster_border + (EXPECTED_DUPLICATION_LENGTH - TSDlength);
				else
					endRightPredictionBorder = left_mate_cluster_border;
			} else
				endRightPredictionBorder = getHighestBorder();
		} else if(hasLeftAlignedSplitCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				endRightPredictionBorder = left_aligned_split_border;
			} else
				endRightPredictionBorder = getHighestBorder();
		}

		//For two mate clusters with duplication,
		//the border is the mate cluster border increased up to half of expected target site duplication length
		//Otherwise the highest coordinate is taken
		else if(hasLeftMateCluster() && hasRightMateCluster()){
			if(hasTSD().equals(DUPLICATION)){
				int TSDlength = getTSDlength();
				if(TSDlength < EXPECTED_DUPLICATION_LENGTH)
					endRightPredictionBorder = left_mate_cluster_border + (EXPECTED_DUPLICATION_LENGTH - TSDlength) / 2;
				else
					endRightPredictionBorder = left_mate_cluster_border;
			} else
				endRightPredictionBorder = getHighestBorder();
		}

		//For single clusters add the fuzziness
		//unless the size of the right mate cluster is unexpectedly small. Then increase it to the expected size
		else if(hasLeftAlignedSplitCluster()){
			endRightPredictionBorder = this.left_aligned_split_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasRightAlignedSplitCluster()){
			endRightPredictionBorder = this.right_aligned_split_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}else if(hasRightMateCluster()){
			endRightPredictionBorder = this.right_mate_cluster_border + SINGLE_CLUSTER_BORDER_FUZZINESS;
		}
		else{
			if(this.left_cluster_length >= this.max_expected_cluster_size - SINGLE_CLUSTER_BORDER_FUZZINESS){
				endRightPredictionBorder = (this.left_mate_cluster_border + SINGLE_CLUSTER_BORDER_FUZZINESS);
			}
			else
				endRightPredictionBorder = (this.left_mate_cluster_border + this.max_expected_cluster_size - this.left_cluster_length);
		}

		return endRightPredictionBorder;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MobilePrediction other = (MobilePrediction) obj;
		if (discordant_read_names == null) {
			if (other.discordant_read_names != null)
				return false;
		} else if (!discordant_read_names.equals(other.discordant_read_names))
			return false;
		if (left_aligned_split_border != other.left_aligned_split_border)
			return false;
		if (left_aligned_split_hits != other.left_aligned_split_hits)
			return false;
		if (left_cluster_length != other.left_cluster_length)
			return false;
		if (left_mate_cluster_border != other.left_mate_cluster_border)
			return false;
		if (left_mate_hits != other.left_mate_hits)
			return false;
		if (mobile_mappings == null) {
			if (other.mobile_mappings != null)
				return false;
		} else if (!mobile_mappings.equals(other.mobile_mappings))
			return false;
		if (original_reference == null) {
			if (other.original_reference != null)
				return false;
		} else if (!original_reference.equals(other.original_reference))
			return false;
		if (right_aligned_split_border != other.right_aligned_split_border)
			return false;
		if (right_aligned_split_hits != other.right_aligned_split_hits)
			return false;
		if (right_cluster_length != other.right_cluster_length)
			return false;
		if (right_mate_cluster_border != other.right_mate_cluster_border)
			return false;
		if (right_mate_hits != other.right_mate_hits)
			return false;
		if (sample_names == null) {
			if (other.sample_names != null)
				return false;
		} else if (!sample_names.equals(other.sample_names))
			return false;
		if (split_read_names == null) {
			if (other.split_read_names != null)
				return false;
		} else if (!split_read_names.equals(other.split_read_names))
			return false;
		return true;
	}

	public void determineSomatic(Properties props) {

		//Check whether a normal predictions file has been provided, then the output is automatically somatic
		if (props.containsKey(MobileDefinitions.NORMAL_PREDICTIONS)
				&& FileValidation.fileValid(props.getProperty(MobileDefinitions.NORMAL_PREDICTIONS))) {
			this.output_somatic = true;
			this.somatic = SOMATIC;
			return;
		}

		//Check whether this is valid normal/tumor calling
		ArrayList<String> allSamples = new ArrayList<String>(Arrays.asList(props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP,0)));
		if(allSamples.size() != 2) return;
		if(!props.containsKey(MobileDefinitions.MINIMUM_NORMAL_SUPPORTING_READS)) return;
		if(props.containsKey(MobileDefinitions.NORMAL_TUMOR_CALLING) && "true".equals(props.getProperty(MobileDefinitions.NORMAL_TUMOR_CALLING).toLowerCase()))
			this.output_somatic = true;
		else
			return;

		//Classify the prediction as somatic/germline based on the presence of a single sample
		String normalSampleName = props.getProperty(MobileDefinitions.NORMAL_SAMPLE);
		ArrayList<String> predictionSamples = new ArrayList<String>(sample_names);
		int minNormalReads = Integer.parseInt(props.getProperty(MobileDefinitions.MINIMUM_NORMAL_SUPPORTING_READS));
		if(predictionSamples.size() == 1){
			if(predictionSamples.get(0).equals(normalSampleName))
				somatic = GERMLINE;
			else
				somatic = SOMATIC;
		}
		//When it has reads from both samples, classify it based on the number of reads
		else {
			int normalSampleCount = sample_counts.get(normalSampleName);
			if(normalSampleCount >= minNormalReads)
				somatic = GERMLINE;
			else if(normalSampleCount > 0)
				somatic = ARTIFACT;
			else
				somatic = SOMATIC;
		}
	}

	public String getSomatic() {
		return somatic;
	}

	public int[] getLeftMateClusterBounderies(){
		if(left_cluster_sam_record == null)
			return null;

		int[] bounderies = new int[2];
		bounderies[0] = left_cluster_sam_record.getAlignmentStart()-1;
		bounderies[1] = left_cluster_sam_record.getAlignmentEnd();
		return bounderies;
	}
	public int[] getRightMateClusterBounderies(){
		if(right_cluster_sam_record == null)
			return null;

		int[] bounderies = new int[2];
		bounderies[0] = right_cluster_sam_record.getAlignmentStart()-1;
		bounderies[1] = right_cluster_sam_record.getAlignmentEnd();
		return bounderies;
	}
	public int[] getSplitClusterBounderies(){
		if(samrecord_splitcluster == null)
			return null;

		int[] bounderies = new int[2];

		if(hasLeftAlignedSplitCluster()){
			String splitEndsString = samrecord_splitcluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_RIGHTCLIPPED_ENDS).toString();
			int[] splitEnds = Arrays.stream(splitEndsString
					.replace("[", "").replace("]","")
					.split(",",0)).mapToInt(Integer::parseInt).toArray();
			for(int splitEnd: splitEnds){
				if(bounderies[0] == 0 || splitEnd < bounderies[0])
					bounderies[0] = splitEnd;
				if(bounderies[1] == 0 || splitEnd > bounderies[1])
					bounderies[1] = splitEnd;
			}
		}
		if(hasRightAlignedSplitCluster()){
			String splitEndsString = samrecord_splitcluster.getAttribute(MobileDefinitions.SAM_TAG_SPLIT_LEFTCLIPPED_ENDS).toString();
			int[] splitEnds = Arrays.stream(splitEndsString
					.replace("[", "").replace("]","")
					.split(",",0)).mapToInt(Integer::parseInt).toArray();
			for(int splitEnd: splitEnds){
				if(bounderies[0] == 0 || splitEnd < bounderies[0])
					bounderies[0] = splitEnd;
				if(bounderies[1] == 0 || splitEnd > bounderies[1])
					bounderies[1] = splitEnd;
			}
		}

		return bounderies;
	}
}