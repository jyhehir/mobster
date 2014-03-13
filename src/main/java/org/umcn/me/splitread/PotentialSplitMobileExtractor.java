package org.umcn.me.splitread;

import net.sf.samtools.*;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.apache.commons.cli.*;
import org.umcn.gen.sam.QualityProcessing;
import org.umcn.gen.sam.SAMSilentReader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This class contains methods to extract reads which are hard-trimmed
 * when mapped to the reference sequence. The trimmed (clipped) parts
 * of these reads are outputted as double-encoded fastq files (bwa) or csfasta
 * and qual files (mosaik) [extract_option == w]. Alternatively the whole sequence
 * of the reads which align with hard-trimming can be outputted to
 * double-encoded fastq files or csfasta and qual files [extract_option == c].
 * Only those hard-clipped reads are outputted which have a minimum
 * number of bases trimmed (specified by the user and by the min_hard_clipping
 * variable in this class) 
 * 
 * This class is designed specifically for colorspace single-end reads mapped with
 * either lifescope or bioscope which allows for hardtrimmed mapping. This
 * class should however also work for colorspace single-end reads which are
 * mapped with other software which also allows for hard-trimming.
 * 
 * TODO 
 * - move makeDoubleEncodeMap, doubleEncode and decodeQualString to org.umcn.me.util package
 * - check whether a bam file gets outputted
 * - when mosaik is chosen as mapping_program min_nr_high_qual_bases & min_quality
 *   are not taken into account
 * - both fastq and (csfasta and qual) are outputted now [either fastq or (csfasta and qual) will
 *   stay empty though]. Only one or the other should be outputted based on the chosen mapping program
 * - make use of the ClippedRead class which was developed after this class
 *  
 * @author Djie
 *
 */
public class PotentialSplitMobileExtractor {
	
	
	private static int min_hard_clipping = 20; //min number of hard clipping on 1 side for reads to
	                                          //be extracted
	
	private static Logger logger = Logger.getLogger("PotentialSplitMobileExtractor");
	
	private static final String[] extract_options = {"c", "w"}; //extract Whole read or clipped sequence
																//for reads which meet the min_hard_clipping
																//requirement
	
	private static final String[] mapping_programs = {"bwa", "mosaik"}; //program which will be used for the
																	    //mobile mapping part
	
	private static int min_nr_high_qual_bases = 14; //a hard-trimmed read with min_hard_clipping
													//will only be extracted when it has this nr of min
													//high quality bases
	private static int min_quality = 10; //high quality base has a phred score >= this variable

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String input = "";
		String prefix = "";
		String extractAlgorithm = "";
		String mappingProgram = "";
		Options options;
		
		BasicConfigurator.configure();
		
		HelpFormatter formatter = new HelpFormatter();
		options = createOptions();	
		
		if (args.length == 0){
			formatter.printHelp("java -Xmx4g -jar PotentialSplitMobileExtractor.jar [-OPTIONS]" , options);

		}else {
			CommandLineParser parser = new GnuParser();
			try {
				CommandLine line = parser.parse(options, args);
				input = line.getOptionValue("in");
				prefix = line.getOptionValue("out");
				min_hard_clipping = Integer.parseInt(line.getOptionValue("minh", Integer.toString(min_hard_clipping)));
				
				mappingProgram = line.getOptionValue("tool", mapping_programs[0]);
				extractAlgorithm = line.getOptionValue("extr", extract_options[0]);
				
				logger.info("This program will assume Sanger (Phred+33) quality encoding!");
				logger.info("Extracting reads with minimal hard clipping of: " + min_hard_clipping);
				logger.info("Extracted reads algorithm: " + extractAlgorithm);
				logger.info("Will prepare output for: " + mappingProgram);
				logger.info("Only extracting reads where the clipped end has at least " +
						min_nr_high_qual_bases + " bases with phred of " + min_quality);
				
				extractPotentialSplitReads(input, prefix, extractAlgorithm, mappingProgram);
				
			} catch (ParseException e) {
				logger.fatal("Error in parsing CLI arguments: " + e.getMessage());
			} catch (IOException e){
				logger.fatal("IO error in extractPotentialSplitReads: " + e.getMessage());
				logger.fatal(e.getStackTrace());
			}
		}
		
		BasicConfigurator.configure();
	}
	
	/**
	 * Create the cmd line options for this class
	 * @return
	 */
	private static Options createOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("BAM file containing hard-clipped reads");
		
		options.addOption(OptionBuilder.create("in"));
		
		OptionBuilder.withArgName("int");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Minimum amount of hard clipping to regard " +
		"read as potential split-mobile read, default: " + min_hard_clipping);
		
		options.addOption(OptionBuilder.create("minh"));
				
		OptionBuilder.withArgName("prefix");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output prefix for csfasta & qual or fastq files containing potential mobile split-reads");
		
		options.addOption(OptionBuilder.create("out"));
		
		OptionBuilder.withArgName("char");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Extract whole read sequence (w) or clipped read sequence (c) " +
				"default is: " + extract_options[0]);
		
		options.addOption(OptionBuilder.create("extr"));
		
		OptionBuilder.withArgName("string");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Which alignment program will be used for mapping against mobile ref." +
				"Options are mosaik or bwa. Default is: " + mapping_programs[0]);
		
		options.addOption(OptionBuilder.create("tool"));
		
		return options;
	}
	
	/**
	 * @return returns a map containing color codes as keys
	 * and their corresponding double-encoded bases as values
	 */
	public static Map<Character, String> makeDoubleEncodeMap(){
		
		Map<Character, String> doubleEncode = new HashMap<Character, String>();
		doubleEncode.put('0', "A");
		doubleEncode.put('1', "C");
		doubleEncode.put('2', "G");
		doubleEncode.put('3', "T");
		doubleEncode.put('.', "N");
				
		return doubleEncode;
	}
	
	/**
	 * Double encode a colorspace sequence. BWA for instance requires
	 * double encoded sequences for colorspace aligning. 
	 * 
	 * @param encodeMap a map containing color codes as keys and their
	 * corresponding double-encoded bases as values
	 * @param seq a colorspace sequence to be encoded
	 * @return a string containing the double-encoded equivalent of the colorspace sequence
	 */
	public static String doubleEncode(Map<Character, String> encodeMap, String seq){
		
		StringBuilder newSeq = new StringBuilder();
		
		for(int i = 0; i < seq.length(); i++){
			newSeq.append(encodeMap.get(seq.charAt(i)));
		}
		
		return newSeq.toString();
	}
	
	/**
	 * This function extracts reads which are clipped for a minimum nr of bases
	 * as specified by min_nr_high_qual_bases on one side. Only the clipped
	 * sequence or the whole sequence is extracted when the read matches
	 * the minimum nr of bases threshold and the min_nr_high_qual_bases &
	 * min_quality threshold (note this is not yet implemented for MOSAIK).
	 * Extracted reads are outputted to fastq and csfasta qual formats
	 * @param in original BAM file location with hard-trimmed reads
	 * @param out output prefix for csfasta, qual, fastq and bam files
	 * @param extractMethod extract only clipped part of sequence (c)
	 * or the whole read (w)
	 * @param tool bwa or mosaik, which tool to use for mobile mapping
	 * @throws IOException
	 */
	public static void extractPotentialSplitReads(String in, String out,
			String extractMethod, String tool) throws IOException{
		File inBam = new File(in);
		//File outBam = new File(out + "_potentialsplitmob.bam");
		PrintWriter outFa = new PrintWriter(new FileWriter(out + "_potentialsplitmobF3.csfasta"), true);
		PrintWriter outQual = new PrintWriter(new FileWriter(out + "_potentialsplitmobF3_QV.qual"), true);
		PrintWriter outFq = new PrintWriter(new FileWriter(out + "_potentialsplitmob.fastq"), true);
		Matcher match = null;
		Pattern begHardClipPattern;
		String cigar;
		String hardClipCigar;
		String hardClipCigar2;
		StringBuilder header = new StringBuilder();
		StringBuilder pattern = new StringBuilder();
		String readSeq = "";
		String readQual = "";
		int hardClipNr;
		int hardClipNr2;
		int nrHits;
		int c = 0;
		int d = 0;
		boolean begClip; //whether clipping occurs at beginning or end
		
		Map<Character, String> encodeMap = makeDoubleEncodeMap();
		
		SAMFileReader inputSam = new SAMSilentReader(inBam);
		//SAMFileHeader headerSam = inputSam.getFileHeader();
		//SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerSam,
				//false, outBam);
		
		Pattern hardClipPattern = Pattern.compile("\\d+H"); //at least 1 digit followed by H, which is for Hard Clipping

		
		for(SAMRecord rec : inputSam ){
			cigar = rec.getCigarString();
			nrHits = Integer.parseInt(rec.getAttribute("NH").toString());
			//Only get uniquely aligned reads
			if (nrHits == 1){
				
				match = hardClipPattern.matcher(cigar);
				
				while (match.find()){
					
					hardClipCigar = match.group();
					
					//hardClipNr == nr of bases which are clipped
					hardClipNr = Integer.parseInt(hardClipCigar.substring(0, hardClipCigar.length() - 1));
					
					if (hardClipNr >= min_hard_clipping){
						header.append(">");
						header.append(rec.getReadName());
						
						//Extract whole read sequence
						if(extractMethod.equals("w")){
							
							if(tool.equals("mosaik")){
								readSeq = rec.getAttribute("CS").toString();
								readQual = decodeQualString(rec.getAttribute("CQ").toString(), 33);
							}else if(tool.equals("bwa")){
								//substring from 2, because you don't want primer base and contaminated
								//adjacent color (bwa and bowtie mentiont this specifically
								readSeq = doubleEncode(encodeMap, rec.getAttribute("CS").toString().substring(2));
								readQual = rec.getAttribute("CQ").toString();
							}
							
							break;
							//extract only clipped sequence when extractMethod == c
						}else if(extractMethod.equals("c")){
							//Check whether read is hard-clipped on two sides; use the side with largest hard clip
							if (match.find()){
								hardClipCigar2 = match.group();
								hardClipNr2 = Integer.parseInt(hardClipCigar2.substring(0, hardClipCigar2.length() - 1));
								if (hardClipNr2 > hardClipNr){
									begClip = false;
									hardClipNr = hardClipNr2;
								}else{
									begClip = true;
								}
							//read is hard clipped on one side, determine on which side:
							}else{
								pattern.append("^");
								pattern.append(hardClipCigar);
								begHardClipPattern = Pattern.compile(pattern.toString());
								if (begHardClipPattern.matcher(cigar).find()){
									begClip = true;
								}else{

									begClip = false;
								}
								pattern.setLength(0);
							}
							boolean negStrand = rec.getReadNegativeStrandFlag();
							
							//in these cases the hard-clipped part is the beginning of the original read seq
							if ((negStrand && !begClip) || (!negStrand && begClip)){
								
								if (tool.equals("mosaik")){
									readSeq = rec.getAttribute("CS").toString().substring(0, hardClipNr + 1); //+1 (incl adaptor base)
									readQual = decodeQualString(rec.getAttribute("CQ").toString().substring(0, hardClipNr), 33);
								}else if(tool.equals("bwa")){
									//Trim primer base + first color
									readSeq = rec.getAttribute("CS").toString().substring(2, hardClipNr + 1);
									readSeq = doubleEncode(encodeMap, readSeq);
									readQual = rec.getAttribute("CQ").toString().substring(1, hardClipNr);
								}

							//in these cases the hard-clipped part is the end of the original read seq
							}else if((negStrand && begClip) || (!negStrand && !begClip)){
								readSeq = rec.getAttribute("CS").toString();
								readSeq = readSeq.substring(readSeq.length() - hardClipNr, readSeq.length());
								readQual = rec.getAttribute("CQ").toString();
								
								if (tool.equals("mosaik")){
									//+1, because MOSAIK has to have qual length - 1 of seq length 
									readQual = readQual.substring(readQual.length() - hardClipNr + 1, readQual.length());
									readQual = decodeQualString(readQual, 33);
								}else if(tool.equals("bwa")){
									readQual = readQual.substring(readQual.length() - hardClipNr, readQual.length());
									readSeq = doubleEncode(encodeMap, readSeq);
								}
								

							}
						}
						c++;
						
						if (tool.equals("mosaik")){
							outFa.println(header.toString());
							outQual.println(header.toString());
							outFa.println(readSeq);
							outQual.println(readQual);
						}else if (tool.equals("bwa")){
							if(QualityProcessing.passesQualityCheck(QualityProcessing.decodeQualStringToVector(readQual, 33),
									min_nr_high_qual_bases, min_quality)){
								outFq.println(header.toString());
								outFq.println(readSeq);
								outFq.println("+");
								outFq.println(readQual);
								d++;
							}
						}
						//outputSam.addAlignment(rec);
						header.setLength(0);
					}
				}
			}
		}

		logger.info("Nr of reads with clipped end: " + c);
		logger.info("Nr or reads that pass Qual check: " + d);
		inputSam.close();
		//outputSam.close();
		outFa.close();
		outQual.close();
		outFq.close();
				
	}
	
	/**
	 * Decode an encoded qual String like ddf`feed]`]_Ba_^__[YBBB
	 * to a numeric qual String as specified in the .qual format
	 * @param eQual: the encoded qual string
	 * @param zero: the ASCII number which equals to a quality of 0.
	 * For instance 33 for Sanger encoding.
	 * @return Numeric qual string.
	 */
	public static String decodeQualString(String eQual, int zero){
		StringBuilder dQual = new StringBuilder();
		
		for (int i = 0; i < eQual.length(); i++){
			dQual.append((int) eQual.charAt(i) - zero);
			dQual.append(" ");
		}
		
		return dQual.toString().trim();
	}
}
