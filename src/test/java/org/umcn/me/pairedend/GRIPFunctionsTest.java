package org.umcn.me.pairedend;

import java.io.File;
import java.util.Vector;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.BasicConfigurator;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.SAMSilentReader;

import junit.framework.TestCase;

public class GRIPFunctionsTest extends TestCase {

	
	private Vector<MobilePrediction> predictions;
	private File grips = new File(this.getClass().getResource("/test_data/A105_testcase_GRIP_clusters.sam").getFile());
	
	
	static{
		BasicConfigurator.configure();
	}
	
	/**
	 * Construct prediction vector with:
	 *  1x SET
	 *  2x SETBP1
	 *  4x ZKSCAN3
	 *  8x GBP4
	 *  3x SKA3
	 */
	public void setUp(){
		this.predictions = this.readInRecords();
		
	}
	
	public void testReducingPredictionsWithMultipleSourceGenes(){
		
		Vector<MobilePrediction> filteredPredictions = GRIPFunctions.reducePredictionsBasedOnSource(this.predictions, 3);
		
		assertEquals(3, filteredPredictions.size());
		
		filteredPredictions = GRIPFunctions.reducePredictionsBasedOnSource(this.predictions, 4);
		
		assertEquals(6, filteredPredictions.size());
		
		filteredPredictions = GRIPFunctions.reducePredictionsBasedOnSource(this.predictions, 5);
		
		assertEquals(10, filteredPredictions.size());
		
	}
	
	private Vector<MobilePrediction> readInRecords(){
		
		Vector<MobilePrediction> preds = new Vector<MobilePrediction>();
		SAMSilentReader reader = new SAMSilentReader(this.grips);
		
		for (SAMRecord rec : reader){
			try {
				preds.add(new MobilePrediction(100, 100, 100, rec));
			} catch (IllegalSAMPairException e) {
				e.printStackTrace();
				fail();
			}
		}
		
		reader.close();
		return preds;
	}
	
	
}
