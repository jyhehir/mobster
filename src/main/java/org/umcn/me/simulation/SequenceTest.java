package org.umcn.me.simulation;

import java.util.*;

import org.apache.log4j.BasicConfigurator;
import org.umcn.gen.sequence.InvalidSequenceException;
import org.umcn.gen.sequence.Sequence;

import junit.framework.TestCase;

/**
 * Testclass for both MobileSequence.java and TargetSequence.java
 * 
 * NOTE: MAKE SURE ASSERTIONS ARE ENABLED FOR VM
 * No AssertionErrors should be thrown.
 * @author Djie
 *
 */
public class SequenceTest extends TestCase{


	public static void main(String[] args) {
		BasicConfigurator.configure();
		//Do not change underneath statements otherwise some assertions will fail
		String testSeq = "AAAATTTT";
		String testSeq2 = "AGGCTAGGCATGACCT";
		List<int[]> testcoordinates = new ArrayList<int[]>();
		List<int[]> testexclusions = new ArrayList<int[]>();
		List<int[]> toolongexclusions = new ArrayList<int[]>();
		List<int[]> wrongorderexcl = new ArrayList<int[]>();
		Map<Integer, Sequence> insertMap = new HashMap<Integer, Sequence>();
		Map<Integer, Sequence> wrongMap = new HashMap<Integer, Sequence>();
		Map<Integer, Sequence> wrongMap2 = new HashMap<Integer, Sequence>();
		Map<Integer, Sequence> wrongMap3 = new HashMap<Integer, Sequence>();
		
		int[] x = {1,5};
		int[] y = {10, 15};
		testcoordinates.add(x);
		testcoordinates.add(y);
		
		int[] excl1 = {4, 10};
		int[] excl2 = {16, 20};
		int[] toolongexcl1 = {4, 6, 8};
		int[] wrongorder = {6, 4};
		testexclusions.add(excl1);
		testexclusions.add(excl2);
		toolongexclusions.add(toolongexcl1);
		toolongexclusions.add(excl2);
		wrongorderexcl.add(wrongorder);
		
		try {
			MobileSequence testMob = new MobileSequence(testSeq);
			
			//MobileSequence tests
			//--------------------
			
			//Test getters and setters of MobileSeq object
			testMob.setName("L1HS");
			testMob.setFamilyName("L1");
			testMob.setDescription("Active LINE-1 element");
			assert testMob.getName().equals("L1HS");
			assert testMob.getDescription().equals("Active LINE-1 element");
			assert testMob.getFamilyName().equals("L1");
			
			//Test truncate method
			assert testMob.truncate(false, 2); // Truncate last 2 nucleotides
			assert testMob.getSequence().equals("AAAATT");
			assert !testMob.truncate(true, 20); //should fail as truncation > seq.length			
			assert testMob.truncate(true, 2); //truncate first 2 nucleotides
			assert testMob.getSequence().equals("AATT"); //This should be resulting seq after 2bp trimmer from both ends
			assert !testMob.truncate(true, 4); //This should fail as truncation of whole seq is not allowed
			assert testMob.getSequence().equals("AATT"); //Seq should therefore not have been changed
			assert testMob.getLeftTruncationSize() == 2; 
			assert testMob.getRightTruncationSize() == 2;
			assert testMob.truncate(true, 1); //1bp is still availabe from left
			assert testMob.getSequence().equals("ATT");
			assert testMob.getLeftTruncationSize() == 3; //Truncationsize should get updated
			assert testMob.truncate(false, 1); //1bp is still available from right
			assert testMob.getRightTruncationSize() == 3; //Truncationsize should get updated
		
			//Test inversion methods and subsequent isValidInsertion method
			MobileSequence testMob2 = new MobileSequence(testSeq2);
			assert testMob2.inverseSubSequence(0,4);
			assert testMob2.getSequence().equals("CGGATAGGCATGACCT");
			assert !testMob2.inverseSubSequence(0,1); //Should fail as inversion of 1nt does not make sense
			assert !testMob2.inverseSubSequence(2, 6); //Should fail as first inversion coordinate overlaps with previous inversion event
			assert testMob2.inverseSubSequence(8, 12); //This should work
			assert testMob2.getSequence().equals("CGGATAGGGTACACCT");
			assert !testMob2.inverseSubSequence(5, 9); //Should fail end coordinate overlaps with previous inversion event
			assert !testMob2.inverseSubSequence(11, 30); //Should fail end coordinate > seq.length
			assert !testMob2.inverseSubSequence(-1, 13); //Should fail begin coordinate < 0
			assert testMob2.getInversionCoordinates().size() == 2;
			System.out.println(Arrays.toString(testMob2.getInversionCoordinates().get(0)));
			System.out.println(Arrays.toString(testMob2.getInversionCoordinates().get(1)));
			assert testMob2.inverseSubSequence(4, 7); //Test inserting inversion event in middle of two insertion events
			assert testMob2.getInversionCoordinates().size() == 3; //Now 3 inversion events have taken place.
			
			//End of MobileSequence tests
			//---------------------------
			
			

			//TargetSequence tests
			//-----------------------
			TargetSequence<MobileSequence> testTar = new TargetSequence<MobileSequence>("CGATGGCCATAC");
			testTar.setName("chr20");
			assert testTar.getName().equals("chr20"); //assert set and getName is correct
			
			//Test getTSD method without exclusion regions
			assert testTar.getTSD(10, 2).getSequence().equals("AT");
			assert testTar.getTSD(4, 10).getSequence().equals("CGAT");
			assert testTar.getTSD(3, -2) == null;
			assert testTar.getTSD(-3, -3) == null;
			assert testTar.getTSD(15, 3) == null;
			assert testTar.getTSD(-2, 3) == null;
			
			//Test isValidInsertion
			assert testTar.insertSingle(0, new Sequence("TTTT")); //Inserting sequences in TS without exclusion regions but with window of 100 should still be legal
			assert testTar.getSequence().equals("TTTTCGATGGCCATAC");
			
			TargetSequence<MobileSequence> testTar2 = new TargetSequence<MobileSequence>("CGGGATGGCCCAAAGTAAACCAGT"); //seq.length 24
			testTar2.setExclusionRegions(testexclusions); //4-10, 16-20
			testTar2.setExclusionWindow(0); //set to 0 for now as we are testing with short seq
			int[] testje = {13, 15};
			int[] testje2 = {23, 24};
			testTar2.addExclusionRegion(testje);
			assert testTar2.addExclusionRegion(testje2); //setting an exclusion region with highest coordinate == seq.length is acceptable
			
			assert testTar2.getExclusionWindow() == 0;
			assert testTar2.isValidInsertion(0); //Inserting at beginning: concatenating two sequences is legal
			assert !testTar2.isValidInsertion(5); //Illegal as exclusion region is 4-10
			assert !testTar2.isValidInsertion(4); //Illegal as exclusion region is 4-10
			assert testTar2.isValidInsertion(24);//Inserting at end is legal (concatenating two sequences)
			assert testTar2.isValidInsertion(2); //Legal
			assert testTar2.isValidInsertion(15); //Legal
			assert testTar2.isValidInsertion(20); //Legal because exclusion region is 16-20 (20 is non-inclusive)
			
			//Test insertion methods with and without exclusion regions
			assert testTar2.insertSingle(0, new Sequence("TTTT")); //Insertion is legal
			assert testTar2.getSequence().equals("TTTTCGGGATGGCCCAAAGTAAACCAGT"); //Thus sequence should be this
			assert testTar2.getOriginalSequenceLength() == 24 ; //But original length should not be updated
			assert testTar2.insertSingle(24, new Sequence("NNNN")) == false; //Second insertion is not (yet) legal (use insertMultiple)
			assert testTar2.getSequence().equals("TTTTCGGGATGGCCCAAAGTAAACCAGT");
			
			TargetSequence<Sequence> testTar3 = new TargetSequence<Sequence>("CGGGATGGCCCAAAGTAAACCAGT");
			assert !testTar3.setExclusionRegions(toolongexclusions); //int[] array of size 3 should not be accepted
			assert !testTar3.setExclusionRegions(wrongorderexcl); //first number in int[] may not be higher than 2nd number.
			assert testTar3.getExclusionRegions().size() == 0; //Failed setting of exclusion regions should result in empty exclusion list
			assert testTar3.setExclusionRegions(testexclusions); //4-10, 16-20
			assert testTar3.getExclusionRegions().size() == 2;
			insertMap.put(0, new Sequence("NNNN"));
			insertMap.put(3, new Sequence("NNNN"));
			insertMap.put(22, new Sequence("NNNNN"));
			assert testTar3.setExclusionWindow(0); //Default is 100, set to 0 for short test sequence
			assert testTar3.insertMultiple(insertMap, true); //should be able to set as all insertion coordinates are legit
			assert !testTar3.isInsertable(); //after multiple insertions, TS should not be insertable anymore
			assert testTar3.getSequence().equals("NNNNCGGNNNNGATGGCCCAAAGTAAACCANNNNNGT");
			assert !testTar3.insertMultiple(insertMap, true); //Multiple insert can not occur anymore
			assert !testTar3.insertSingle(12, new Sequence("TTTT"));
			
			TargetSequence<Sequence> testTar4 = new TargetSequence<Sequence>("CGGGATGGCCCAAAGTAAACCAGT");
			testTar4.setExclusionWindow(0);
			wrongMap.put(0, new Sequence("NNNN"));
			wrongMap.put(55, new Sequence("NNNN"));
			assert !testTar4.insertMultiple(wrongMap, true); //as 55 is longer than seq.length this should fail
			assert testTar4.isInsertable(); //After failed insertion, this sequence should still be insertable
			wrongMap2.put(0, new Sequence("NNNN"));
			wrongMap2.put(-4, new Sequence("NNNN"));
			assert !testTar4.insertMultiple(wrongMap2, true); //should fail negative index not allowed
			assert testTar4.setExclusionRegions(testexclusions); //4-10, 16-20
			wrongMap3.put(1, new Sequence("NNNN"));
			wrongMap3.put(6, new Sequence("NNNN"));
			assert !testTar4.insertMultiple(wrongMap3, true); //should fail, 6 is in exclusion region
			assert testTar4.insertMultiple(insertMap, true); //should work (0, 3, 20)
			
			//MobileSequence copy tests
			MobileSequence testMob4 = new MobileSequence(testMob);
			assert testMob4.getFamilyName().equals("L1");
			testMob.setFamilyName("Alu");
			assert testMob4.getFamilyName().equals("L1");
			assert testMob.getFamilyName().equals("Alu");


			
			
		}catch (InvalidSequenceException ex) {
			ex.getMessage();
		}

	}

}
