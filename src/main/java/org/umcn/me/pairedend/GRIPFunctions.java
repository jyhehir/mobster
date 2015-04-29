package org.umcn.me.pairedend;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.apache.log4j.Logger;

public class GRIPFunctions {
	
	public static Logger logger = Logger.getLogger(GRIPFunctions.class);

	/**
	 * Only return predictions which occur a maximum of x times based on the name of the source genes in the vector.
	 * @param predictions
	 * @param max an integer: only predictions with < max same source genes are kept
	 * @return a copy of predictions containing only the predictions which adhere to the filtering criteria.
	 */
	public static Vector<MobilePrediction> reducePredictionsBasedOnSource(Vector<MobilePrediction> predictions, int max){
		
		Map<String, Integer> sourceGeneCounts = new HashMap<String, Integer>();
		Set<String> sourceGenesToLog = new HashSet<String>();
		Vector<MobilePrediction> predictionsToKeep = new Vector<MobilePrediction>();
		sourceGeneCounts = getSourceGeneCounts(predictions);
		
	
		for (int i=0; i < predictions.size(); i++){
			MobilePrediction pred = predictions.get(i);
			
			for (String sourceGene : pred.getMobileMappings()){
				if (sourceGeneCounts.containsKey(sourceGene)){
					int actualCount = sourceGeneCounts.get(sourceGene);
					
					if (actualCount < max){
						predictionsToKeep.add(pred);
					}					
					if (sourceGenesToLog.add(sourceGene)){
						logger.debug("Source gene: " + sourceGene +  " - occurs: " + actualCount);
					}
				}
			}
		}
		
		logger.info("Separate predictions are allowed to have a maximum of: " + max + " same source genes");
		logger.info("Number of predictions filtered because of same source gene: " + (predictions.size() - predictionsToKeep.size()));
		
		return predictionsToKeep;


	}
	
	public static Map<String, Integer> getSourceGeneCounts(Vector<MobilePrediction> predictions){
		
		Map<String, Integer> counts = new HashMap<String, Integer>();
		
		for (MobilePrediction pred : predictions){
			for (String sourceGene : pred.getMobileMappings()){
				
				if (counts.containsKey(sourceGene)){
					int newCount = counts.get(sourceGene) + 1;
					counts.put(sourceGene, newCount);
				}else{
					counts.put(sourceGene, 1);
				}
				
			}
		}
		return counts;
		
	}
	
}
