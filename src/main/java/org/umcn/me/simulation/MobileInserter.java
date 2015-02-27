package org.umcn.me.simulation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileWriter;

import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.RichSequence;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.bio.seq.io.FastaHeader;

import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.apache.commons.cli.*;
import org.umcn.gen.sequence.FastaSequence;
import org.umcn.gen.sequence.InvalidSequenceException;
import org.umcn.gen.sequence.Sequence;
import org.umcn.me.parsers.FastaParserFast;

/**
 * Class for simulating Mobile Element Insertion events in DNA Sequences (TargetSequences)
 * DNA sequences can be inserted in inclusion regions only (based on .bed file for instance
 * only in exons) or "whole-genomically".
 * Insertions with characteristics (like size of TSD) are written to output file
 * 
 * TODO
 * rewrite of writeMobileInsertions (use StringBuilder for efficiency)
 * implement that sequences (MobileSequences) can truncate or invert when inserted
 * 
 * @author Djie
 *
 */
public class MobileInserter {

	public static Logger logger = Logger.getLogger("MobileInserter");

	private static final int rep_mask_begin = 6; // index of beginning genomic
													// coordinate of element in
													// repmask file
	private static final int rep_mask_end = 7; // index of end genomic
												// coordinate of element in
												// repmask file

	private static String[] me_families = { "ALU", "SVA", "L1", "HERV" };


	public static void runMobileInserter(String targetFa, String excl, String incl, String insertFa,
			MobileInserterProperties props, String outFa, String outFile) throws IOException,
			InvalidSequenceException {
		List<int[]> exclRegions = new ArrayList<int[]>();
		Vector<FastaSequence> targetSeqVec = null; //target sequence vector
		Vector<FastaSequence> meSeqVec = null; // mobile elements vector
		TargetSequence<MobileSequence> targetSeq = null; //Target Sequence
		
		//Vector of Target Sequences with insertions
		Vector<FastaSequence> targetWithInsertsVec = new Vector<FastaSequence>();
		Map<String, List<MobileSequence>> meMap = new HashMap<String, List<MobileSequence>>();
		Map<Integer, MobileSequence> insertMap;
		List<Integer> randomInserts;
		List<Integer> inclusionRegions = null;

		boolean insertSucces = false;
		
		//Print out Mobile Inserter Properties used
		logger.info(props.toString());

		targetSeqVec = readFastaBioJava(targetFa);
		meSeqVec = readFastaBioJava(insertFa);
		
		if (targetSeqVec == null || targetSeqVec.size() != 1) {
			logger.error("Target fasta may only contain 1 sequence as of yet. Terminating.");
			return;
		}
		
		targetSeq = new TargetSequence<MobileSequence>(targetSeqVec.get(0).getSequence());

		// TODO: call function to check whether MobileInserterProperties
		// settings are correct
		targetSeq.setName(targetSeqVec.get(0).getId().trim());
		targetSeq.setExclusionWindow(props.getExclusionWindow());

		if (excl == "" || excl == null) {
			logger.info("No exclusion file given, no exclusions will be set for target sequence");
		} else {
			exclRegions = getExclusionsFromRepMask(excl);
			if (exclRegions == null) {
				logger.error("Exclusion file provided, but no exclusions found "
						+ "check format. Terminating");
				return;
			}
			targetSeq.setExclusionRegions(exclRegions);
			logger.info("Nr of excl regions: " + targetSeq.getExclusionRegions().size());
		}

		if (incl == "" || incl == null){
			logger.info("No inclusion regions set.");
			randomInserts = createRandomInsertions(targetSeq, props);
		}else{
			inclusionRegions = new ArrayList<Integer>(getInclusionFromBed(incl, props.getInclusionWindow()));
			randomInserts = createRandomInsertions(targetSeq, props, inclusionRegions);
		}

		meMap = createMobileSequences(meSeqVec, props);
		logger.info("Nr of excl regions after random insertion" + targetSeq.getExclusionRegions().size());
		// System.out.println(randomInserts);
		insertMap = createMobileInsertMap(randomInserts, meMap, targetSeq, props);
		// System.out.println(insertMap);
		insertSucces = targetSeq.insertMultiple(insertMap, false);
		logger.info("insertSucces" + insertSucces);
		// System.out.println(targetSeq);
		if (insertSucces) {
			//targetSeq.setName(targetSeq.getName() + "_insertions");
			writeMobileInsertions(insertMap, targetSeq.getName(), outFile);
			targetWithInsertsVec.add(new FastaSequence(targetSeq.getName(), targetSeq.getSequence()));
			FastaParserFast.writeFasta2(targetWithInsertsVec, outFa);
		}

	}

	public static void main(String[] args) {
		Options options;
		String exclFile = "";
		String targetFile = "";
		String insertFile = "";
		String outFile = "";
		String outFasta = "";
		String inclFile = null;
		BasicConfigurator.configure();
		options = makeCmdOptions();
		CommandLineParser parser = new GnuParser();
		
		try {
			CommandLine line = parser.parse(options, args);
			if (line.hasOption("incl")){
				inclFile = line.getOptionValue("incl");
			}
			if (line.hasOption("excl")){
				exclFile = line.getOptionValue("excl");
			}
			targetFile = line.getOptionValue("target");
			insertFile = line.getOptionValue("ins");
			outFile = line.getOptionValue("outt");
			outFasta = line.getOptionValue("outf");
			MobileInserterProperties properties = new MobileInserterProperties();
			runMobileInserter(targetFile, exclFile, inclFile, insertFile, properties, outFasta, outFile);
		} catch (ParseException e1) {
			logger.error("Error parsing CLI arguments: " + e1.getMessage());
		} catch (IOException e) {
			logger.error("Error in running runMobileInserter " + e.getMessage());
		} catch (InvalidSequenceException e) {
			logger.error("Error in running runMobileInserter " + e.getMessage());
		}

	}
	
	private static Options makeCmdOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName(".out");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("RepMask exclusion file");
		
		options.addOption(OptionBuilder.create("excl"));
		
		OptionBuilder.withArgName("fa");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Reference fasta file");
		
		options.addOption(OptionBuilder.create("target"));
		
		OptionBuilder.withArgName("fa");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("fasta file containing sequences to insert");
		
		options.addOption(OptionBuilder.create("ins"));
		
		OptionBuilder.withArgName("fa");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("output fasta file");
		
		options.addOption(OptionBuilder.create("outf"));
		
		OptionBuilder.withArgName("txt");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("output text file");
		
		options.addOption(OptionBuilder.create("outt"));
		
		OptionBuilder.withArgName("txt");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("inclusion file");
		
		options.addOption(OptionBuilder.create("incl"));
		
		return options;
	}

	public static Map<String, List<MobileSequence>> addTsdsToMobileSeqs(
			Map<String, List<MobileSequence>> meMap, MobileInserterProperties props) {
		List<MobileSequence> mobileList = new ArrayList<MobileSequence>();
		int tsdNr = 0;

		for (String family : meMap.keySet()) {
			mobileList = meMap.get(family);
			Collections.shuffle(mobileList); // pick MEs randomly for TSD

			if (family.equals("ALU")) {
				tsdNr = props.getAluTsdNrs();
			} else if (family.equals("L1")) {
				tsdNr = props.getL1TsdNrs();
			} else if (family.equals("SVA")) {
				tsdNr = props.getSVATsdNrs();
			} else if (family.equals("HERV")) {
				tsdNr = props.getHervTsdNrs();
			}
			for (int i = 0; i < tsdNr; i++) {
				MobileSequence me = mobileList.get(i);
				me.setInsertWithTsd(true);
				mobileList.set(i, me);
			}
			meMap.put(family, mobileList);
		}
		return meMap;
	}

	public static Map<Integer, Sequence> createInsertMap(List<Integer> insertList,
			Map<String, List<? extends Sequence>> meMap) {
		Map<Integer, Sequence> insertMap = new HashMap<Integer, Sequence>();
		int counter = 0;

		for (String family : meMap.keySet()) {
			for (Sequence seq : meMap.get(family)) {
				insertMap.put(insertList.get(counter), seq);
				counter++;
			}
		}
		if (insertList.size() != insertMap.keySet().size()) {
			logger.warn("Size of insertion positions list not same as number of insertion seqs");
			return null;
		}
		return insertMap;
	}

	public static Map<Integer, MobileSequence> createMobileInsertMap(List<Integer> insertList,
			Map<String, List<MobileSequence>> meMap, TargetSequence<MobileSequence> target,
			MobileInserterProperties props) {
		Map<Integer, MobileSequence> insertMap = new HashMap<Integer, MobileSequence>();
		int counter = 0;
		int tsdSize;
		Random ran = new Random();
		int minTsd = props.getMinTsdSize();
		int maxTsd = props.getMaxTsdSize();
		int insertion;

		for (String family : meMap.keySet()) {
			for (MobileSequence seq : meMap.get(family)) {
				insertion = insertList.get(counter);
				if (seq.insertsWithTsd()) {
					tsdSize = ran.nextInt(maxTsd - minTsd + 1) + minTsd;
					Sequence tsd = target.getTSD(insertion, tsdSize);
					seq.setTsdSize(tsd.getLength());
					seq.insert(seq.getSequence().length(), tsd);
				}
				insertMap.put(insertion, seq);
				counter++;
			}
		}
		if (insertList.size() != insertMap.keySet().size()) {
			logger.warn("Size of insertion positions list not same as number of ME seqs");
			return null;
		}
		return insertMap;
	}

	/**
	 * Function which categorizes Mobile Sequences into families (i.e. ALU, L1)
	 * and modifies them according to the MobileInserterProperties supplied:
	 * - randomly expand/decrease the amount of mobile elements to the nr
	 * supplied in the properties file
	 * - add TSDs to specified amount of mobile elements
	 * - still TODO add Inversions, truncations
	 * @param seqs
	 * @param props
	 * @return
	 */
	public static Map<String, List<MobileSequence>> createMobileSequences(
			Vector<FastaSequence> seqs, MobileInserterProperties props) {
		Map<String, List<MobileSequence>> meSeqs = new HashMap<String, List<MobileSequence>>();
		for (FastaSequence seq : seqs) {
			try {
				MobileSequence me = new MobileSequence(seq.getSequence());
				for (String family : me_families) {
					if (seq.getHeader().toUpperCase().indexOf(family) != -1) {
						if (meSeqs.containsKey(family)) {
							meSeqs.get(family).add(me);
						} else {
							List<MobileSequence> meList = new ArrayList<MobileSequence>();
							meList.add(me);
							meSeqs.put(family, meList);
						}
						me.setFamilyName(family);
						me.setName(seq.getId());
					}
				}
			} catch (InvalidSequenceException e) {
				// should not occur
				logger.warn(e.getMessage());
			}
		}
		for (String family : me_families){
			if(meSeqs.get(family) != null){
				logger.info("Nr of " + family +
						" elements found: " + meSeqs.get(family).size());
			}else{
				logger.info("0 elements found for " + family);
			}
		}
		// Move this section to a check mobileProps function
		if (props.getAluInsertionNrs() > 0 && meSeqs.get("ALU") == null) {
			System.out
					.println("0 ALU sequences found but user specified ALU sequences should be inserted");
			return null;
		} else if (props.getL1InsertionNrs() > 0 && meSeqs.get("L1") == null) {
			System.out
					.println("0 L1 sequences found but user specified L1 sequences should be inserted");
			return null;
		} else if (props.getSVAInsertionNrs() > 0 && meSeqs.get("SVA") == null) {
			System.out
					.println("0 SVA sequences found but user specified SVA sequences should be inserted");
			return null;
		} else if (props.getHervTsdNrs() > 0 && meSeqs.get("HERV") == null){
			System.out
			.println("0 HERV sequences found but user specified SVA sequences should be inserted");
		}

		//increase or decrease number of me's per family according to props file
		meSeqs = createSpecifiedNrMobileSequences(meSeqs, props);
		if (props.getTotalTsds() > 0){
			logger.info("Adding TSDs to Mobile Sequences");
			meSeqs = addTsdsToMobileSeqs(meSeqs, props);
		}
		
		//add here future methods for mobile element modifications (truncations, inversions)
		
		
		return meSeqs;
	}

	public static List<Integer> createRandomInsertions(TargetSequence<? extends Sequence> target,
			MobileInserterProperties properties) {
		Random ran = new Random();
		List<Integer> insertionList = new ArrayList<Integer>();
		int[] exclRegion = new int[2];
		int totalInsertions = properties.getTotalInsertions();

		logger.info("Creating random insertions without inclusion regions specified.");
		for (int i = 0; i < totalInsertions; i++) {
			boolean succes = false;
			while (!succes) {
				int insertion = ran.nextInt(target.getOriginalSequenceLength());
				if (target.isValidInsertion(insertion)) {
					insertionList.add(insertion);
					exclRegion[0] = insertion;
					exclRegion[1] = insertion + 1;

					target.addExclusionRegion(exclRegion);
					succes = true;
				}
				
			}
		}
		return insertionList;
	}
	
	public static List<Integer> createRandomInsertions(TargetSequence<? extends Sequence> target,
			MobileInserterProperties properties, List<Integer> inclusions) {
		Random ran = new Random();
		List<Integer> insertionList = new ArrayList<Integer>();
		int[] exclRegion = new int[2];
		int totalInsertions = properties.getTotalInsertions();

		logger.info("Creating random insertions with inclusion regions specified.");
		
		for (int i = 0; i < totalInsertions; i++) {
			boolean succes = false;
			while (!succes) {
				int insertionIndex = ran.nextInt(inclusions.size());
				int insertion = inclusions.get(insertionIndex);
				if (target.isValidInsertion(insertion)) {
					insertionList.add(insertion);
					exclRegion[0] = insertion;
					exclRegion[1] = insertion + 1;
					target.addExclusionRegion(exclRegion);
					succes = true;
				}
				
			}
		}
		return insertionList;
	}

	private static Map<String, List<MobileSequence>> createSpecifiedNrMobileSequences(
			Map<String, List<MobileSequence>> meMap, MobileInserterProperties props) {
		
		int nrElements = 0;
		List<MobileSequence> meList = null;
		Random ran = new Random();

		for (String family : meMap.keySet()) {
			List<MobileSequence> tempList = new ArrayList<MobileSequence>();
			meList = meMap.get(family);
			if (family.equals("ALU")) {
				nrElements = props.getAluInsertionNrs();
			} else if (family.equals("SVA")) {
				nrElements = props.getSVAInsertionNrs();
			} else if (family.equals("L1")) {
				nrElements = props.getL1InsertionNrs();
			} else if (family.equals("HERV")){
				nrElements = props.getHervInsertionNrs();
			}

			for (int i = 0; i < nrElements; i++) {
				int randomInt = ran.nextInt(meList.size());
				try {
					MobileSequence me = new MobileSequence(meList.get(randomInt));
					tempList.add(me);
					//System.out.println(me.getFamilyName());
				} catch (InvalidSequenceException e) {
					// should not occur
					logger.warn("Not a valid sequence in createSpecifiedNrMobileSequence");
				}

			}
			meMap.put(family, tempList);

		}
		return meMap;
	}

	public static List<int[]> getExclusionsFromRepMask(String repfile) throws IOException {
		
		String line = "";
		String[] lineSplit;
		int begin;
		int end;
		List<int[]> exclusions = new ArrayList<int[]>();
		boolean header = true;

		BufferedReader br = new BufferedReader(new FileReader(repfile));

		while ((line = br.readLine()) != null) {
			if (!header) {
				try {
					lineSplit = line.split("\t");
					begin = Integer.parseInt(lineSplit[rep_mask_begin]);
					end = Integer.parseInt(lineSplit[rep_mask_end]);
					int[] exclusion = { begin, end };
					exclusions.add(exclusion);
				} catch (Exception e) {
					// Catches both NumberFormatException when casting
					logger.warn("Error in parsing repeatmasker file for getting exclusion regions "
							+ "; probably due to wrong file format or wrong specified columns. "
							+ e.getMessage());
					br.close();
					return null;
				}
			}else{
				header = false;
			}
		}
		br.close();
		return exclusions;
	}

	//TODO: build StringBuilder into this method for efficiency
	public static void writeMobileInsertions(Map<Integer, MobileSequence> meMap, String target,
			String out) {
		PrintWriter outFile;
		String sep = "\t";
		
		Integer[] insertions = meMap.keySet().toArray(new Integer[0]);
		Arrays.sort(insertions);

		try {
			outFile = new PrintWriter(new FileWriter(out, true), true); // appends
																		// ==
																		// true

			// write header
			outFile.println("target" + sep + "insertPosition" + sep + "meName" + sep + "meFamily"
					+ sep + "insertionLength" + sep + "tsd" + sep + "tsdSize" + sep + "inversion");

			for (int insertion : insertions) {
				MobileSequence me = meMap.get(insertion);
				outFile.print(target + sep + Integer.toString(insertion) + sep + me.getName() + sep
						+ me.getFamilyName() + sep
						+ Integer.toString(me.getLength() - me.getTsdSize()) + sep
						+ Boolean.toString(me.insertsWithTsd()) + sep + me.getTsdSize() + sep);
				for (int[] inversion : me.getInversionCoordinates()) {
					outFile.print(Integer.toString(inversion[0]) + "-");
					outFile.print(Integer.toString(inversion[1]) + ",");
				}
				outFile.println();
			}
			outFile.close();
		} catch (IOException e) {
			logger.error(e.getMessage());
		}
	}

	public static void writeSeqWithInsertions(Vector<TargetSequence<MobileSequence>> seqVec,
			String out) {
		String header;
		String seq;
		FileOutputStream fos;
		
		File file = new File(out);
		FastaHeader fasHeader = new FastaHeader();
		
		file.delete();
		//Do not write version and namespace info to fasta header
		fasHeader.setShowVersion(false);
		fasHeader.setShowNamespace(false);
		
		try {
			fos = new FileOutputStream(file, true);

			for (TargetSequence<MobileSequence> s : seqVec) {
				header = s.getName();
				seq = s.getSequence();
				try {
					RichSequence.IOTools.writeFasta(fos, DNATools.createDNASequence(seq, header),null,  fasHeader);
					
				} catch (IllegalSymbolException e){
					
				} catch (IOException e){
					logger.warn("IOException when writing seqs to fasta: " + e.getMessage());
				}
			}
			fos.close();
		} catch (IOException e) {
			logger.warn("File not found when writing seqs to fasta: " + e.getMessage());
		}
	}

	// maybe move to FastaParser
	public static Vector<FastaSequence> readFastaBioJava(String fa) throws IOException {
		RichSequence rec;
		Vector<FastaSequence> seqs = new Vector<FastaSequence>();
		BufferedReader br = new BufferedReader(new FileReader(fa));
		RichSequenceIterator iterator = RichSequence.IOTools.readFastaDNA(br, null);
		while (iterator.hasNext()){
			try {
				rec = iterator.nextRichSequence();
				//may be replaced by StringBuilder method for efficiency
				seqs.add(new FastaSequence(">" + rec.getName(), rec.seqString()));
			} catch (InvalidSequenceException e) {
				logger.warn("Invalid sequence found" + " in "  + fa + ", " + e.getMessage());
			} catch (BioException e){
				logger.warn("Invalid sequence found" + " in "  + fa + ", " + e.getMessage());
			}
		}
		logger.info("Read " + String.valueOf(seqs.size()) + " sequences from " + fa);
		return seqs;
	}
	
	public static Set<Integer> getInclusionRegions(String infile, int exclusion) throws IOException{
		BufferedReader br;
		String line;
		boolean header = true;
		int exonStarts = 9;
		int exonEnds = 10;
		String[] starts;
		String[] ends;
		Set<Integer> inclusionSet = new HashSet<Integer>();
		
			br = new BufferedReader(new FileReader(infile));
			while ((line = br.readLine()) != null ){
				if (!header){
					starts = line.split("\t")[exonStarts].replaceAll(",$", "").split(","); //remove trailing comma
					ends = line.split("\t")[exonEnds].replaceAll(",$", "").split(",");
					for (int i = 0; i < starts.length; i++){
						int start = Integer.parseInt(starts[i]) + exclusion;
						int end = Integer.parseInt(ends[i]) - exclusion ;
						if (start < end){
							for (int j = start; j < end; j ++){
								inclusionSet.add(j);
							}
						}
					}
				}else{
					header = false;
				}			
			}
		br.close();
		logger.info("Nr of nucleotides that are insertable: " +inclusionSet.size());
		return inclusionSet; 
	}
	
	public static Set<Integer> getInclusionFromBed(String infile, int exclusion) throws IOException{
		BufferedReader br;
		String line;

		Set<Integer> inclusionSet = new HashSet<Integer>();
		
			br = new BufferedReader(new FileReader(infile));
			while ((line = br.readLine()) != null ){
				String[] split = line.split("\t", -1);
				int start = Integer.parseInt(split[1]) + exclusion;
				int end = Integer.parseInt(split[2]) - exclusion;
				if (start < end) {
					for (int i = start; i < end; i++) {
						inclusionSet.add(i);
					}
				}
				
			}
			
			br.close();
			logger.info("Nr of nucleotides that are insertable: " + inclusionSet.size());
			return inclusionSet;
			
	}
}
