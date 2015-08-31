package org.umcn.me.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.log4j.Logger;
import org.umcn.gen.region.InvalidRegionException;
import org.umcn.gen.region.LabeledRegion;
import org.umcn.me.pairedend.MobilePrediction;

public class CollectionUtil {

	public static Logger logger = Logger.getLogger("CollectionUtil");
	
	public static <T extends Object> String toString(Collection<T> collection){
		StringBuilder sb = new StringBuilder();
		
		if(collection.size() == 0){
			return "";
		}
		
		for (T object : collection){
			sb.append(object.toString());
			sb.append(",");
		}
		
		return sb.substring(0, sb.length() - 1);
	}
	
	public static <T extends Object> String toString(Collection<T> collection, String delim){
		StringBuilder sb = new StringBuilder();
		
		if(collection.size() == 0){
			return "";
		}
		
		for (T object : collection){
			sb.append(object.toString());
			sb.append(delim);
		}
		
		return sb.substring(0, sb.length() - 1);
	}
	
	public static Map<LabeledRegion, MobilePrediction> mobilePredictionsToRegions(Vector<MobilePrediction> preds,
			int insecurity){
		
		Map<LabeledRegion, MobilePrediction> map  = new HashMap<LabeledRegion, MobilePrediction>();
		
		int c = 0;
		
		for (MobilePrediction pred : preds){
			c++;
			int start = pred.getLeftPredictionBorder() + 1 - insecurity;
			int end = pred.getRightPredictionBorder() + insecurity;
			
			if (start < 1){
				start = 1;
			}
			
			String chr = pred.getChromosome();
			String label = Integer.toString(c);
			try {
				map.put(new LabeledRegion(label, chr, start, end), pred);
			} catch (InvalidRegionException e) {
				logger.error("Could not convert MobilePrediction to LabeledRegion!");
			}
		}
		
		return map;
		
	}
	
	/**
	 * Convert a collection of readnames to a map with readname strings as key and readname object as values.
	 * Note: if there are duplicate readnames, only the last readname object with the same readname string
	 * will be added to the map.
	 * @param readnames
	 * @return
	 */
	public static Map<String, ReadName> readNamesToMap(Collection<ReadName> readnames){
		
		Map<String, ReadName> map = new HashMap<String, ReadName>();
		
		for (ReadName name : readnames){
			map.put(name.readName, name);
		}
		return map;	
	}
	
	/**
	 * Helper function to keep track how many times a string occurs.
	 * 
	 * If key occurs in map then add +1 to the corresponding integer value of key.
	 * If key does not occur in map add integer value of 1 to corresponding key.
	 * @param key
	 * @param countMap
	 */
	public static void addKeyToCountMap(String key, Map<String, Integer> countMap){
		
		if (countMap.containsKey(key)){
			countMap.put(key, countMap.get(key) + 1);
		}else{
			countMap.put(key, 1);
		}
	}
	
	public static String mapToString(Map<? extends Object,? extends Object> map){
		String mapString = map.toString();
		return mapString.substring(1, mapString.length() - 1);
	}
	
	public static Set<String> returnOverlappingKeysInMap(Map<String, ? extends Object> map1, Map<String, ? extends Object> map2){
		
		Set<String> map1Keys = map1.keySet();
		Set<String> map2Keys = map2.keySet();
		Set<String> mapCopy = new HashSet<String>(map1Keys);
		mapCopy.retainAll(map2Keys);

		return mapCopy;
	}
	
}
