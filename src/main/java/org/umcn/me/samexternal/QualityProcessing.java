package org.umcn.me.samexternal;

import java.math.BigDecimal;
import java.util.Vector;


/**
 * Collection of short methods which mostly deal
 * with QUAL string and PHRED operations. 
 * 
 * @author Djie
 *
 */
public class QualityProcessing {

	
	/**
	 * Decode an encoded qual String like ddf`feed]`]_Ba_^__[YBBB
	 * to a numeric qual String as specified in the .qual format
	 * @param eQual: the encoded qual string
	 * @param zero: the ASCII number which equals to a quality of 0.
	 * For instance 33 for Sanger / Illumina 1.8+ encoding.
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
	
	/**
	 * Decode an encoded qual String like ddf`feed]`]_Ba_^__[YBBB
	 * to a numeric qual String as specified in the .qual format
	 * @param eQual: the encoded qual string
	 * @param zero: the ASCII number which equals to a quality of 0.
	 * For instance 33 for Sanger / Illumina 1.8+ encoding.
	 * @return Numeric qual vector.
	 */
	public static Vector<Integer> decodeQualStringToVector(String eQual, int zero){
		Vector<Integer> dQual = new Vector<Integer>();
		
		for (int i = 0; i < eQual.length(); i++){
			dQual.add((int) eQual.charAt(i) - zero);
		}
		
		return dQual;
	}
	
	/**
	 * @param qualityVector vector with phred scores (numeric)
	 * @param nrBases nr of bases that have to have phred score higher or equal than param quality
	 * @param quality phred score
	 * @return
	 */
	public static boolean passesQualityCheck(Vector<Integer> qualityVector,
			int nrBases, int quality){
		
		int sum = 0;
		
		for (int i = 0; i < qualityVector.size(); i++){
			if (qualityVector.get(i) >= quality ){
				sum += 1;
			}
		}
		
		if(sum >= nrBases){
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * This method calculates the median of a SORTED integer vector
	 * @param estimates vector containing integers
	 * @return the median
	 */
	public static int calculateMedian(Vector<Integer> estimates) {
		int length;

		length = estimates.size();
		
		if (length % 2 == 0){
			
			int number1 = estimates.get((length / 2) - 1);
			int number2 = estimates.get((length / 2));
			return ((number1 + number2) /  2);
		}else{
			return estimates.get((length - 1) / 2);
		}

	}
	
	/**
	 * Calculate an average over a vector of integers with 2 digits
	 * precision
	 * @param integers vector of integers to calculate average upon
	 * @return average with 2 digits of precision
	 */
	public static double calculateAverage(Vector<Integer> integers){
		int sum = 0;
		double length = integers.size();
		double average = 0.0;
		
		for(int i : integers){
			sum += i;
		}
		
		average = sum / length;
		
		return (double) Math.round(average * 100) / 100;
		
	}
	
	/**
	 * Round a double
	 * @param unrounded double to be rounded
	 * @param precision number of digits precision
	 * @param roundingMode usually set to BigDecimal.ROUND_HALF_UP
	 * @return
	 */
	public static double round(double unrounded, int precision, int roundingMode)
	{
	    BigDecimal bd = new BigDecimal(unrounded);
	    BigDecimal rounded = bd.setScale(precision, roundingMode);
	    return rounded.doubleValue();
	}

	/**
	 * Creates a String DNA sequence of custom length
	 * containing only base N.
	 * @param number length of sequence to be returned.
	 * @return
	 */
	public static String createNSequence(int number){
		StringBuilder nseq = new StringBuilder();
		for(int i=0; i < number; i++){
			nseq.append("N");
		}
		return nseq.toString();
	}


}
