//TODO: still needs to be converted to be compatible with BioInfoBase2 and new NGS package

//package org.umcn.me.splitread;
//
//import java.io.File;
//import java.io.FileWriter;
//import java.io.IOException;
//import java.io.PrintWriter;
//import java.util.Collections;
//import java.util.HashSet;
//import java.util.Set;
//import java.util.Vector;
////import java.util.regex.Matcher;
////import java.util.regex.Pattern;
//
//import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMRecord;
//import net.sf.samtools.SAMRecordIterator;
//
//import org.apache.commons.cli.CommandLine;
//import org.apache.commons.cli.CommandLineParser;
//import org.apache.commons.cli.GnuParser;
//import org.apache.commons.cli.HelpFormatter;
//import org.apache.commons.cli.OptionBuilder;
//import org.apache.commons.cli.Options;
//import org.apache.commons.cli.ParseException;
//import org.apache.log4j.Logger;
//import org.apache.log4j.BasicConfigurator;
//import org.umcn.maf.feature.IntFeature;
//import org.umcn.maf.feature.StringFeature;
//import org.umcn.maf.region.AnnotatedRegion;
//import org.umcn.maf.region.AnnotatedRegionInterface;
//import org.umcn.maf.region.AnnotatedRegionSet;
//import org.umcn.maf.region.InvalidRegionException;
//import org.umcn.maf.region.RegionComparator;
//import org.umcn.maf.region.RegionInterface;
//import org.umcn.gen.sam.SAMSilentReader;
//import org.umcn.me.util.*;
//import org.umcn.ngs.annotaters.ExonIntronAnnotater;
//import org.umcn.me.sam.InvalidCategoryException;
//import org.umcn.me.splitread.ClippedRead;
//
////NOTE ALL PREDICTIONS ARE 1-BASED INSTEAD OF UCSC 0-BASED AND THE PAIRED-END PREDICTIONS
///**
// * 
// * This class contains methods to cluster anchors from a split-read method
// * together and add all kinds of annotation to these clusters (==predictions),
// * upon which these clusters can later be filtered, such as average phred base
// * quality of clipped ends.
// * 
// * 
// * Cluster all anchors together which are within x bp together (max_distance).
// * x is calculated from the alignmentstart:
// * So for a x of 15
// * 
// * xxxxxxxxxxxxxxx
// * -------| (first clipped anchor)
// * ATCGATTAGCCCCATTTGAGACT
// * 
// * All reads which fall even partially in the xxxxxx region will also be clustered.
// * Regardless whether they map to the same mobile element family or not. This can
// * not be done accurately with short length reads.
// * 
// * After clustering predictions (AnnotatedRegion instance) are made containing
// * all sorts of annotation like average phred score of clipped end,
// * distance between left clipped ends etc. are added.
// * 
// * The prediction window for a prediction ranges from
// * Leftmost clipped position minus prediction_border
// * Rightmost clipped position plus prediction_border
// * 
// * The estimate for insertion position is the median of all clipped positions
// * 
// * All predictions which overlap with a repeat are filtered. [in this case all repeats,
// * because we use exome data and clipped short ends can easily be mistaken for
// * different sorts of repeats]
// * 
// * @author Djie
// *
// */
//public class SplitAnchorClusterer {
//
//	public static Logger logger = Logger.getLogger("SplitAnchorClusterer");
//	
//	private static int max_distance = 50; //maximum distance between anchor MEI estimates
//	private static int min_support = 1; //nr of minimum supporting reads
//	private static int prediction_border = 15; //this number of bases is added to both sides of estimate
//	private static int min_length_support = 20;
//	private static String repmask_location = "./hg19_reps_correctorder.out";
//	private static String refgene_location = "./refGene_hg19.txt";
//	
//	//private final static String DB_RESOURCE_FILE = "./lib/ucsc_db_hg19.properties";
//	private final static String COLUMN_NRHITS = "nr_hits";
//	private final static String COLUMN_MOBILE = "mobile_category";
//	private final static String COLUMN_ESTIMATE = "estimate";
//	private final static String COLUMN_READS = "reads";
//	private final static String COLUMN_SAMPLE = "sample";
//	private final static String COLUMN_SIDES = "both_sides";
//	private final static String COLUMN_LENGTH = "avg_mob_length";
//	private final static String COLUMN_LENGTH_LEFT = "avg_mob_length_left";
//	private final static String COLUMN_LENGTH_RIGHT = "avg_mob_length_right";
//	private final static String COLUMN_AVG_QUAL = "avg_clipped_qual";
//	private final static String COLUMN_AVG_QUAL_LEFT = "avg_clipped_qual_left";
//	private final static String COLUMN_AVG_QUAL_RIGHT = "avg_clipped_qual_right";
//	private final static String COLUMN_FRAC_IDEN_LEFT = "fraction_sameclipping_left";
//	private final static String COLUMN_FRAC_IDEN_RIGHT = "fraction_sameclipping_right";
//	private final static String COLUMN_NR_LEFT = "nr_leftclipped";
//	private final static String COLUMN_NR_RIGHT = "nr_rightclipped";
//	private final static String COLUMN_MAXDIST_LEFT = "max_clipdistance_left";
//	private final static String COLUMN_MAXDIST_RIGHT = "max_clipdistance_right";
//	private final static String COLUMN_NR_ALU = "nr_alu";
//	private final static String COLUMN_NR_SVA = "nr_sva";
//	private final static String COLUMN_NR_L1 = "nr_l1";
//
//	
//	
//	public static void main(String[] args) {
//		
//		Options options;
//		String anchorBam;
//		String outPrediction;
//		String index;
//		
//		BasicConfigurator.configure();
//		HelpFormatter formatter = new HelpFormatter();
//		options = createOptions();
//		
//		if (args.length == 0){
//			formatter.printHelp("java -Xmx4g -jar SplitAnchorClusterer.jar", options);
//		}else{
//			CommandLineParser parser = new GnuParser();
//			try {
//				CommandLine line = parser.parse(options, args);
//				
//				anchorBam = line.getOptionValue("in");
//				outPrediction = line.getOptionValue("out") + "_predictions.txt";
//
//				index = line.getOptionValue("index");
//				
//				max_distance = Integer.parseInt(line.getOptionValue("maxdist", Integer.toString(max_distance)));
//				min_support = Integer.parseInt(line.getOptionValue("minsup", Integer.toString(min_support)));
//				min_length_support = Integer.parseInt(line.getOptionValue("minlength", Integer.toString(min_length_support)));
//				repmask_location = line.getOptionValue("rep", repmask_location);
//				refgene_location = line.getOptionValue("refgene", refgene_location);
//				
//				
//				logger.info("Using a maximum distance between MEI supporting ends of anchors of " + max_distance);
//				logger.info("A MEI prediction should be supported by a minimum of " + min_support + " anchors.");
//				logger.info("When a prediction has less than " + min_support + " supporting reads, "+ 
//						"a prediction will still be made, when the average length of reads in the cluster is >= " +
//						min_length_support);
//				
//				runClusterer(anchorBam, index, outPrediction);
//				
//			} catch (ParseException e) {
//				logger.fatal("Error in parsing CLI arguments");
//			} catch (IOException e){
//				logger.fatal("Error in writing reading files: " + e.getMessage());
//			} catch (InvalidRegionException e){
//				logger.fatal("Error in creating prediction. Prediction is invalid region: " + e.getMessage());
//			} catch (InvalidHardClipCigarException e) {
//				logger.fatal("Anchor bam file contains reads without hard-clipping!");
//			}
//			
//		}
//		
//	}
//
//	/**
//	 * The do-work function for this class. Clusters all split anchors, then
//	 * filters all predictions with prediction windows which are within repeats.
//	 * The predictions which are left behind are annotated with the exon_intron
//	 * annotator.
//	 * All predictions are subsequently written to a txt file
//	 * 
//	 * @param inBam (split) anchor bam file (output of ExtractSplitAnchors.java) coordinate sorted
//	 * @param bai associated .bai file of inBam
//	 * @param outString file name for the prediction file to be outputted
//	 * @throws IOException
//	 * @throws InvalidRegionException
//	 * @throws InvalidHardClipCigarException
//	 */
//	private static void runClusterer(String inBam, String bai, String outString) throws IOException,
//	InvalidRegionException, InvalidHardClipCigarException {
//		
//		Vector<AnnotatedRegion> predictions = new Vector<AnnotatedRegion>();
//		Vector<AnnotatedRegion> filteredPredictions = new Vector<AnnotatedRegion>();
//		AnnotatedRegionSet<AnnotatedRegion> filteredPredictionsSet;
////		ExonIntronAnnotater<AnnotatedRegion> eia = new ExonIntronAnnotater<AnnotatedRegion>();
//		
//		//STEP1 Cluster all anchors
//		predictions = clusterAndCreatePredictions(inBam, bai);
//		
//		//STEP2 filter for overlap with repeats
//		filteredPredictions = filterAllPredictionsWithinRepeat(predictions);
//		
//		//STEP3 Annotate with exon-intron annotater
//		filteredPredictionsSet = new AnnotatedRegionSet<AnnotatedRegion>(filteredPredictions);
////		if (!eia.annotateRegions(filteredPredictionsSet, refgene_location, true)){
////			logger.warn("Error in annotating filtered predictions with refGene file");
////		}
//		
//		//STEP4 write predictions to file
//		writePredictions(filteredPredictionsSet, outString);
//	
//	}
//
//	/**
//	 * Write a vector of AnnotatedRegionInterfaces to a file.
//	 * @param filteredPredictions vector of predictions (annotatedregioninterfaces)
//	 * @param outString file location for output prediction file
//	 * @throws IOException
//	 */
//	private static void writePredictions(
//			Vector<? extends AnnotatedRegionInterface> filteredPredictions, String outString) throws IOException {
//		
//		String header = "";
//		
//		PrintWriter outFile = new PrintWriter(new FileWriter(outString), true);
//		
//		if (filteredPredictions.size() != 0 ){
//			header = filteredPredictions.get(0).getHeader();
//		
//			outFile.println(header);
//		
//			for (AnnotatedRegionInterface prediction : filteredPredictions){
//				outFile.println(prediction.toString());
//			}
//		}else{
//			outFile.println("No predictions made.");
//		}
//		
//		outFile.close();
//		
//	}
//
//
//	/**
//	 * Filter all predictions (or any other annotatedregion) which based on start and end postion
//	 * overlap with a repeat.
//	 * @param predictions vector containing annotatedregions, in this case specifically
//	 * containing predictions
//	 * @return vector of annotatedregions which do not overlap with a repeat
//	 * @throws IOException
//	 */
//	private static Vector<AnnotatedRegion> filterAllPredictionsWithinRepeat(
//			Vector<AnnotatedRegion> predictions) throws IOException {
//		
//		Vector<AnnotatedRegion> filteredPredictions = new Vector<AnnotatedRegion>();
//		Vector<RegionInterface> repeats = MobileFilterer.createAnnotatedRepeatVector(repmask_location, false);
//		
//		for (AnnotatedRegion r : predictions){
//			if(RegionComparator.getOverlappingRegionsFast( (RegionInterface) r, repeats).size() == 0){
//				filteredPredictions.add(r);
//			}
//			
//		}
//		
//		logger.info(filteredPredictions.size() + " predictions remain after filtering for all rep mask repeats.");
//		
//		return filteredPredictions;
//	}
//
//
//	/**
//	 * Cluster anchors together when they are within x bp (specified by global
//	 * variable (max_distance) of the alignment start of the 1st anchor.
//	 * All clusters (=predictions) with annotation (fraction_same_clipping, avg_phred_score etc.)
//	 * are collected in a vector and returned
//	 * 
//	 * @param inBam (split) anchor bam file (output of ExtractSplitAnchors.java) coordinate sorted
//	 * @param inBai associated .bai file of inBam
//	 * @return
//	 * @throws IOException
//	 * @throws InvalidRegionException
//	 * @throws InvalidHardClipCigarException
//	 */
//	private static Vector<AnnotatedRegion> clusterAndCreatePredictions(String inBam, String inBai) throws IOException,
//		InvalidRegionException, InvalidHardClipCigarException {
//		
//		int alignStart = 0;
//		int c = 0;
//		String readName = "";
//		String reference = "";
//		String mobileCat = "";
//		AnnotatedRegion prediction = null;
//		Vector<AnnotatedRegion> predictionVectors = new Vector<AnnotatedRegion>();
//		
//		Set<String> readsAlreadyClustered = new HashSet<String>();
//		
//		File inFile = new File(inBam);
//		File index = new File(inBai);
//		
//		SAMFileReader input = new SAMSilentReader(inFile, index);
//		
//		for(SAMRecord rec : input){
//			readName = rec.getReadName();
//			reference = rec.getReferenceName();
//			
//			//TODO make container class for ME attribute tag
//			mobileCat = rec.getAttribute("ME").toString().split(";")[2];
//			
//			if(!readsAlreadyClustered.contains(readName)){
//				
//				alignStart = rec.getAlignmentStart();
//				
//				//do the real clustering and adding annotation with this function call
//				//this function will also return information of the SAM record
//				//in the current loop
//				prediction = getOverlappingAnchorandCreatePrediction(inFile, index,
//						reference, mobileCat, alignStart);
//				
//				if (prediction != null) {
//					//add all reads which are clustered from this prediction
//					//to readsAlreadyClustered, so no more clustering
//					//will be done on these reads in future iterations in this loop
//					for (String read : prediction.getFeatureTable()
//							.get(COLUMN_READS).toString().split(";")) {
//						readsAlreadyClustered.add(read);
//					}
//					prediction.getFeatureTable().remove(COLUMN_READS);
//					predictionVectors.add(prediction);
//					c += Integer.parseInt(prediction.getFeatureTable()
//							.get(COLUMN_NRHITS).toString());
//				}
//
//			}
//		}
//		
//		logger.info(predictionVectors.size() + " unfiltered predictions made from " + c + " reads.");
//		input.close();
//		
//		return predictionVectors;
//
//	}
//
//
//	/**
//	 * Cluster all anchors together which are within x bp together (max_distance).
//	 * x is calculated from the alignmentstart:
//	 * So for a x of 15
//	 * 
//	 * xxxxxxxxxxxxxxx
//	 * -------|
//	 * ATCGATTAGCCCCATTTGAGACT
//	 * 
//	 * All reads which fall even partially in the xxxxxx region will also be clustered.
//	 * Regardless whether they map to the same mobile element family or not. This can
//	 * not be done accurately with short length reads.
//	 * 
//	 * After clustering predictions (AnnotatedRegion instance) are made containing
//	 * all sorts of annotation like average phred score of clipped end,
//	 * distance between left clipped ends etc. are added.
//	 * 
//	 * The prediction window for a prediction ranges from
//	 * Leftmost clipped position minus prediction_border
//	 * Rightmost clipped position plus prediction_border
//	 * 
//	 * The estimate for insertion position is the median of all clipped positions
//	 * 
//	 * @param inFile
//	 * @param index
//	 * @param reference
//	 * @param mobileCat
//	 * @param alignStart
//	 * @return
//	 * @throws IOException
//	 * @throws InvalidRegionException
//	 * @throws InvalidHardClipCigarException
//	 */
//	private static AnnotatedRegion getOverlappingAnchorandCreatePrediction(File inFile, File index,
//			String reference, String mobileCat, int alignStart) throws IOException,
//			InvalidRegionException, InvalidHardClipCigarException {
//		
//		ClippedRead rec = null;
//		int border5 = 0;
//		int border3 = 0;
//		int nrHits = 0;
//		int estimate = 0;
//		
//		int nrAlu = 0;
//		int nrSVA = 0;
//		int nrL1 = 0;
//
//
//		//This will contain all reads (== anchors; == clipped anchors) that
//		//are going to be clustered together
//		ClippedReadSet<ClippedRead> clippedSet = new ClippedReadSet<ClippedRead>();
//
//		StringBuilder reads = new StringBuilder();
//		AnnotatedRegion mobileRegion = null;
//		
//		Vector<Integer> estimates = new Vector<Integer>();
//		
//		SAMFileReader input = new SAMSilentReader(inFile, index);
//		SAMRecordIterator iter = null;
//		String sample = input.getFileHeader().getReadGroups().get(0).getSample();
//		
//		
//		iter = input.queryOverlapping(reference, alignStart, alignStart + max_distance);
//		
//		//Iterate the region from alignment start of first anchor till
//		//alignment start + max_distance. All reads which fall partially
//		//into this region will also be taken into account
//		while (iter.hasNext()){
//			rec = new ClippedRead(iter.next(), false, true);
//			
//			reads.append(rec.getSAMRecord().getReadName());
//			reads.append(";");		
//			clippedSet.add(rec);
//
//		}
//		
//		input.close();
//		
//		//make vector of all clipped positions of anchors
//		estimates = clippedSet.getClippedEnds();
//		Collections.sort(estimates);
//		
//		//Create prediction window
//		border5 = estimates.firstElement() - prediction_border;
//		border3 = estimates.lastElement() + prediction_border;
//		
//		nrHits = clippedSet.size();
//		estimate = calculateMedian(estimates);
//		try {
//			nrAlu = clippedSet.getNrOfAluPredictions();
//			nrSVA = clippedSet.getNrOfSVAPredictions();
//			nrL1 = clippedSet.getNrOfL1Predictions();
//		} catch (InvalidCategoryException e) {
//			logger.error("Error in determining number of Alu, L1, SVA:" + e.getMessage());
//		}
//		
//		//TODO these if-elses do not handle ties yet
//		//not very important as mobile identification is
//		//not very accurate with such short reads
//		if(nrAlu > nrSVA && nrAlu > nrL1){
//			mobileCat = "ALU";
//		}else if(nrSVA > nrAlu && nrSVA > nrL1){
//			mobileCat = "SVA";
//		}else if(nrL1 > nrSVA && nrL1 > nrAlu){
//			mobileCat = "L1";
//		}
//
//		
//		double avgLength = clippedSet.getAverageClippedLength();
//		
//		//only make a prediction when a prediction has a minimum number of support
//		//(supporting reads) or when average length of 1 read is >= min length for
//		//clipped sequence
//		//then add all counts of annotation
//		if (nrHits >= min_support || (int) avgLength >= min_length_support) {
//			mobileRegion = new AnnotatedRegion(reference, border5, border3);
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_MOBILE, mobileCat));
//			mobileRegion.getFeatureTable().addFeature(
//					new IntFeature(COLUMN_NRHITS, nrHits));
//			mobileRegion.getFeatureTable().addFeature(
//					new IntFeature(COLUMN_ESTIMATE, estimate));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_SAMPLE, sample));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_READS, reads.substring(0,
//							reads.length() - 1)));
//			if (clippedSet.getClippedLengthsLeft().size() != 0 && clippedSet.getClippedLengthsRight().size() != 0) {
//				mobileRegion.getFeatureTable().addFeature(
//						new StringFeature(COLUMN_SIDES, "true"));
//			} else {
//				mobileRegion.getFeatureTable().addFeature(
//						new StringFeature(COLUMN_SIDES, "false"));
//			}
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_LENGTH, Double
//							.toString(avgLength)));
//
//
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_LENGTH_LEFT, Double
//							.toString(clippedSet.getAverageClippedLengthLeft())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_LENGTH_RIGHT, Double
//							.toString(clippedSet.getAverageClippedLengthRight())));
//			mobileRegion.getFeatureTable().addFeature(new StringFeature(COLUMN_AVG_QUAL,
//					Double.toString(clippedSet.getAverageBaseQualAll())));
//			
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_AVG_QUAL_LEFT, Double
//							.toString(clippedSet.getAverageBaseQualLeftClipped())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_AVG_QUAL_RIGHT, Double
//							.toString(clippedSet.getAverageBaseQualRightClipped())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_FRAC_IDEN_LEFT, Double
//							.toString(clippedSet.getFractionOfLeftClippedReadsWithSamePos())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_FRAC_IDEN_RIGHT, Double
//							.toString(clippedSet.getFractionOfRightClippedReadsWithSamePos())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_NR_LEFT, Integer
//							.toString(clippedSet.getNumberLeftClippedReads())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_NR_RIGHT, Integer
//							.toString(clippedSet.getNumberRightClippedReads())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_MAXDIST_LEFT, Integer
//							.toString(clippedSet.getMaxDistanceLeftClippedReads())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_MAXDIST_RIGHT, Integer
//							.toString(clippedSet.getMaxDistanceRightClippedReads())));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_NR_ALU, Integer
//							.toString(nrAlu)));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_NR_SVA, Integer
//							.toString(nrSVA)));
//			mobileRegion.getFeatureTable().addFeature(
//					new StringFeature(COLUMN_NR_L1, Integer
//							.toString(nrL1)));			
//			
//			return mobileRegion;
//		}else{
//			return null;
//		}
//			
//	}
//
//
//	/**
//	 * This method calculates the median of a SORTED integer vector
//	 * TODO move this method to org.umcn.me.util
//	 * @param estimates vector containing integers
//	 * @return the median
//	 */
//	public static int calculateMedian(Vector<Integer> estimates) {
//		int length;
//
//		length = estimates.size();
//		
//		if (length % 2 == 0){
//			
//			int number1 = estimates.get((length / 2) - 1);
//			int number2 = estimates.get((length / 2));
//			return ((number1 + number2) /  2);
//		}else{
//			return estimates.get((length - 1) / 2);
//		}
//
//	}
//
//
//	/**
//	 * Create command line options for this class
//	 * @return
//	 */
//	private static Options createOptions() {
//		Options options = new Options();
//		
//		OptionBuilder.withArgName("BAM File");
//		OptionBuilder.hasArg();
//		OptionBuilder.isRequired();
//		OptionBuilder.withDescription("Split anchor bam file, coordinate sorted");
//		
//		options.addOption(OptionBuilder.create("in"));
//		
//		OptionBuilder.withArgName("BAI");
//		OptionBuilder.hasArg();
//		OptionBuilder.isRequired();
//		OptionBuilder.withDescription("Split anchor bai (index) file");
//		
//		options.addOption(OptionBuilder.create("index"));
//		
//		OptionBuilder.withArgName("int");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Maximum distance between anchor hard-trimmed endpoints. " +
//				"Hard-trimmed endpoints within maximum distance will be clustered together. " +
//				"Default is: " + max_distance);
//		
//		options.addOption(OptionBuilder.create("maxdist"));
//		
//		OptionBuilder.withArgName("int");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Minimum required of anchors for 1 MEI event." +
//				"Default is: " + min_support);
//		
//		options.addOption(OptionBuilder.create("minsup"));
//		
//		OptionBuilder.withArgName("prefix");
//		OptionBuilder.hasArg();
//		OptionBuilder.isRequired();
//		OptionBuilder.withDescription("Output prefix for prediction txt file.");
//		
//		options.addOption(OptionBuilder.create("out"));
//		
//		OptionBuilder.withArgName("int");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Min average length for a cluster supporting MEI event, when the " +
//				"number of reads in the cluster is lower than specified bij -minsup.");
//		
//		options.addOption(OptionBuilder.create("minlength"));
//		
//		OptionBuilder.withArgName("File");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Repeat masker outfile containing all repeats on hg19. " +
//				"File should be coordinate sorted. Default location: " + repmask_location);
//		options.addOption(OptionBuilder.create("rep"));
//		
//		OptionBuilder.withArgName("File");
//		OptionBuilder.hasArg();
//		OptionBuilder.withDescription("Refgene hg19 file for exon_inton annotator. " +
//				"Default location: " + refgene_location);
//		options.addOption(OptionBuilder.create("refgene"));
//		
//		return options;
//	}	
//	
//}
