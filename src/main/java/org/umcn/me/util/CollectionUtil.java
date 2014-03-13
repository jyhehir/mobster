package org.umcn.me.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.apache.log4j.Logger;
import org.umcn.gen.annotation.AnnotatedRegionSet;
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
	
}
