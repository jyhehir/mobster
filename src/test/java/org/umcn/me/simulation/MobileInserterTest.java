package org.umcn.me.simulation;

import java.io.*;
import java.util.NoSuchElementException;

import org.biojavax.bio.seq.RichSequence;
import org.biojava.bio.BioException;


import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.gen.sequence.InvalidSequenceException;

/**
 * Test class for MobileInserter.java
 * 
 * In essence each testcase makes sure that the sequence with insertions
 * can be reverse engineered to the original sequence based on the output file
 * with insertions
 * 
 * Make sure -ea (for assert statements) is enabled in VM for this class.
 * 
 * TODO make a real junit class from this class
 *
 * 
 * @author Djie
 *
 */
public class MobileInserterTest {
	public static Logger logger;
	

	public static void main(String[] args) throws BioException {
		BasicConfigurator.configure();
		Logger logger = Logger.getLogger("MobileInserterTest");
		logger.warn("blaat");
		assert reverseEngineerTest1();
		assert reverseEngineerTest2();
		assert reverseEngineerTest3();
		assert reverseEngineerTest4();
//		try {
//			System.out.println(reverseEngineerSingle("D:\\testsequences\\targetsequence.fa", "D:\\testsequences\\outTargetWithActive.fa",
//					"D:\\testsequences\\outTargetWithActive.txt"));
//		} catch (NoSuchElementException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}
	
	/**
	 * Test with no TSD, inversions or truncations set
	 * @return whether tests succeeds
	 */
	public static boolean reverseEngineerTest1(){
		String targetFasta = "D:\\testsequences\\targetsequence.fa";
		String exclFile = "D:\\testsequences\\repbase.out";
		String insertFasta = "D:\\testsequences\\insertsequences.fa";
		String outFasta = "D:\\testsequences\\outTestcase1.fa";
		String outFile = "D:\\testsequences\\outTestcase1stats.txt";
		boolean succes = false;;
		
		File file = new File(outFile);
		file.delete(); //Delete previous instance of outFile as runMobileInserter appends to outFile
		
		MobileInserterProperties props = new MobileInserterProperties();
		props.setAluInsertionNrs(3)
		     .setL1InsertionNrs(3)
		     .setSVAInsertionNrs(3)
		     .setAluTsdNrs(0)
		     .setL1TsdNrs(0)
		     .setSVATsdNrs(0)
		     .setAluInversionNrs(0)
		     .setL1InversionNrs(0)
		     .setSVAInversionNrs(0)
		     .setExclusionWindow(0);
		System.out.println(props.toString());
		try {
			MobileInserter.runMobileInserter(targetFasta, exclFile, null, insertFasta, props, outFasta, outFile);
			succes = reverseEngineerSingle(targetFasta, outFasta, outFile);
		} catch (IOException e) {
			logger.warn(e.getMessage());
		} catch (InvalidSequenceException e) {
			logger.warn(e.getMessage());
		} catch (NoSuchElementException e){
			logger.warn(e.getMessage());
		} catch (BioException e){
			logger.warn(e.getMessage());
		}
		return succes;
	}
	
	/**
	 * Test with different TSDs set, but inversion and truncations off
	 * @return
	 */
	public static boolean reverseEngineerTest2(){
		String targetFasta = "D:\\testsequences\\targetsequence.fa";
		String exclFile = "D:\\testsequences\\repbase.out";
		String insertFasta = "D:\\testsequences\\insertsequences.fa";
		String outFasta = "D:\\testsequences\\outTestcase2.fa";
		String outFile = "D:\\testsequences\\outTestcase2stats.txt";
		boolean succes = false;;
		
		File file = new File(outFile);
		file.delete(); //Delete previous instance of outFile as runMobileInserter appends to outFile
		
		MobileInserterProperties props = new MobileInserterProperties();
		props.setAluInsertionNrs(3)
		     .setL1InsertionNrs(3)
		     .setSVAInsertionNrs(3)
		     .setAluTsdNrs(3)
		     .setL1TsdNrs(2)
		     .setSVATsdNrs(1)
		     .setAluInversionNrs(0)
		     .setL1InversionNrs(0)
		     .setSVAInversionNrs(0)
			 .setExclusionWindow(0);
		System.out.println(props.toString());
		try {
			MobileInserter.runMobileInserter(targetFasta, exclFile, null, insertFasta, props, outFasta, outFile);
			succes = reverseEngineerSingle(targetFasta, outFasta, outFile);
		} catch (IOException e) {
			logger.warn(e.getMessage());
		} catch (InvalidSequenceException e) {
			logger.warn(e.getMessage());
		} catch (NoSuchElementException e){
			logger.warn(e.getMessage());
		} catch (BioException e){
			logger.warn(e.getMessage());
		}
		return succes;
	}
	
	/**
	 * Test with inclusion regions set and exclusion regions. No editing events.
	 * @return
	 */
	public static boolean reverseEngineerTest3(){
		String targetFasta = "D:\\testsequences\\targetsequence.fa";
		String exclFile = null;
		String insertFasta = "D:\\testsequences\\insertsequences.fa";
		String outFasta = "D:\\testsequences\\outTestcase3.fa";
		String outFile = "D:\\testsequences\\outTestcase3stats.txt";
		String inclFile = "D:\\testsequences\\inclusionregions.txt";
		boolean succes = false;
		
		File file = new File(outFile);
		file.delete(); //Delete previous instance of outFile as runMobileInserter appends to outFile
		
		MobileInserterProperties props = new MobileInserterProperties();
		props.setAluInsertionNrs(3)
		     .setL1InsertionNrs(3)
		     .setSVAInsertionNrs(4)
		     .setAluTsdNrs(0)
		     .setL1TsdNrs(0)
		     .setSVATsdNrs(0)
		     .setAluInversionNrs(0)
		     .setL1InversionNrs(0)
		     .setSVAInversionNrs(0)
		     .setInclusionWindow(3)
		     .setExclusionWindow(3);
		System.out.println(props.toString());
		try {
			MobileInserter.runMobileInserter(targetFasta, exclFile, inclFile, insertFasta, props, outFasta, outFile);
			succes = reverseEngineerSingle(targetFasta, outFasta, outFile);
		} catch (IOException e) {
			logger.warn(e.getMessage());
		} catch (InvalidSequenceException e) {
			logger.warn(e.getMessage());
		} catch (NoSuchElementException e){
			logger.warn(e.getMessage());
		} catch (BioException e){
			logger.warn(e.getMessage());
		}
		return succes;
	}
	
	/**
	 * Test with exclusion regions set and exclusion windows set
	 * @return
	 */
	public static boolean reverseEngineerTest4(){
		String targetFasta = "D:\\testsequences\\targetsequence.fa";
		String exclFile = "D:\\testsequences\\repbase.out";
		String insertFasta = "D:\\testsequences\\insertsequences.fa";
		String outFasta = "D:\\testsequences\\outTestcase4.fa";
		String outFile = "D:\\testsequences\\outTestcase4stats.txt";
		boolean succes = false;;
		
		File file = new File(outFile);
		file.delete(); //Delete previous instance of outFile as runMobileInserter appends to outFile
		
		MobileInserterProperties props = new MobileInserterProperties();
		props.setAluInsertionNrs(4)
		     .setL1InsertionNrs(4)
		     .setSVAInsertionNrs(5)
		     .setAluTsdNrs(0)
		     .setL1TsdNrs(0)
		     .setSVATsdNrs(0)
		     .setAluInversionNrs(0)
		     .setL1InversionNrs(0)
		     .setSVAInversionNrs(0)
			 .setExclusionWindow(3);
		System.out.println(props.toString());
		try {
			MobileInserter.runMobileInserter(targetFasta, exclFile, null, insertFasta, props, outFasta, outFile);
			succes = reverseEngineerSingle(targetFasta, outFasta, outFile);
		} catch (IOException e) {
			logger.warn(e.getMessage());
		} catch (InvalidSequenceException e) {
			logger.warn(e.getMessage());
		} catch (NoSuchElementException e){
			logger.warn(e.getMessage());
		} catch (BioException e){
			logger.warn(e.getMessage());
		}
		return succes;
	}
	
	/**
	 * Function to reverse engineer a sequence with insertions to an original sequence
	 * on basis of the output statistics file (ie. file with insertion positions and insertion length).
	 * When you can reverse engineer the sequence with insertions to the original sequence,
	 * you know that the insertion position, insertion length, tsd size reported in the output statics
	 * file are correct. 
	 * @throws BioException 
	 * @throws NoSuchElementException 
	 */
	public static boolean reverseEngineerSingle(String originalFasta, String newFasta, String statFile) throws IOException,
			NoSuchElementException, BioException{
		String originalSeq;
		String newSeq;
		String line;
		String delim = "\t";
		StringBuilder reverseEngineered;
		String[] split;
		
		int index;
		int insertLength;
		int tsdSize;
		int colInsert = 1;
		int colInsertLength = 4;
		int colTsd = 6;
		
		boolean header = true;
		
		BufferedReader br = new BufferedReader(new FileReader(statFile));
		BufferedReader brOriginal = new BufferedReader(new FileReader(originalFasta));
		BufferedReader brNew = new BufferedReader(new FileReader(newFasta));
		
		originalSeq = RichSequence.IOTools
								  .readFastaDNA(brOriginal, null)
								  .nextSequence().seqString();
		newSeq = RichSequence.IOTools
		  					 .readFastaDNA(brNew, null)
		  					 .nextSequence().seqString();
		reverseEngineered = new StringBuilder(newSeq);
		
		while ((line = br.readLine()) != null){
			if (!header){
				split = line.split(delim);
				index = Integer.parseInt(split[colInsert]);
				insertLength = Integer.parseInt(split[colInsertLength]);
				tsdSize = Integer.parseInt(split[colTsd]);
				reverseEngineered.replace(index, index + insertLength + tsdSize, "");
			}else{
				header = false;
			}
		}
		br.close();
		

		return originalSeq.equals(reverseEngineered.toString());
	}
}
