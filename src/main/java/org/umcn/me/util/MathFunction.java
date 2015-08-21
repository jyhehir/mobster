package org.umcn.me.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class MathFunction {

	public static Double getMedianFromDoubles(final List<Double> collection){
		int length;

		length = collection.size();
		ArrayList<Double> colCopy = new ArrayList<Double>(collection);
		
		Collections.sort(colCopy);
		
		if (length == 0){
			return null;
		}else if (length == 1){
			return colCopy.get(0);
		}else if (length % 2 == 0){
			
			double number1 = colCopy.get((length / 2) - 1);
			double number2 = colCopy.get((length / 2));
			return ((number1 + number2) /  2);
		}else{
			return colCopy.get((length - 1) / 2);
		}
	}
	
	//TODO: circumvent duplicate code
	public static Integer getMedianFromIntegers(final List<Integer> collection){
		int length;

		length = collection.size();
		ArrayList<Integer> colCopy = new ArrayList<Integer>(collection);
		
		Collections.sort(colCopy);
		
		if (length == 0){
			return null;
		}else if (length == 1){
			return colCopy.get(0);
		}else if (length % 2 == 0){
			
			int number1 = colCopy.get((length / 2) - 1);
			int number2 = colCopy.get((length / 2));
			return ((number1 + number2) /  2);
		}else{
			return colCopy.get((length - 1) / 2);
		}
	}
}
