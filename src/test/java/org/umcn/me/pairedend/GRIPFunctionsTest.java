package org.umcn.me.pairedend;

import java.io.File;
import java.util.Vector;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.BasicConfigurator;
import org.umcn.me.output.FilterPredictions;
import org.umcn.me.samexternal.IllegalSAMPairException;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.util.SimpleRegion;

import junit.framework.TestCase;

public class GRIPFunctionsTest extends TestCase {

	
	private File grips = new File(this.getClass().getResource("/test_data/A105_testcase_GRIP_clusters.sam").getFile());
	private File gripsWithOverlappingClusters = new File(this.getClass().getResource("/test_data/A105a_overlappingclusters_onchr2.sam").getFile());
	
	static{
		BasicConfigurator.configure();
	}
	
	
	public void testReducingPredictionsWithMultipleSourceGenes(){
		Vector<MobilePrediction> predictions;
		/**
		 * Construct prediction vector with:
		 *  1x SET
		 *  2x SETBP1
		 *  4x ZKSCAN3
		 *  8x GBP4
		 *  3x SKA3
		 */
		predictions = this.readInRecords(this.grips);
		FilterPredictions filter = new FilterPredictions(predictions);

		filter.reducePredictionsBasedOnSource(3);
		
		assertEquals(3, filter.getPredictions().size());
		
		filter.reducePredictionsBasedOnSource(4);
		
		assertEquals(6, filter.getPredictions().size());
		
		filter.reducePredictionsBasedOnSource(5);
		
		assertEquals(10, filter.getPredictions().size());
	}
	
	public void testRemovingOverlappingClusters(){
		
		/**
		 * The prediction set contains 23 predictions of which 12 have non-onverlapping prediction windows.
		 */
		Vector<MobilePrediction> predictions = this.readInRecords(gripsWithOverlappingClusters);
		FilterPredictions filter = new FilterPredictions(predictions);

		filter.removeOverlappingPredictions();
		
		//Assert 12 predictions remain now.
		assertEquals(12, filter.getPredictions().size());
		
		/**
		 * 
		 */
		Vector<SimpleRegion> srs = new Vector<SimpleRegion>();
		
		for (MobilePrediction pred : filter.getPredictions()){
			SimpleRegion sr = pred.predictionWindowToRegion();
			srs.add(sr);
		}
		assertEquals(12, srs.size());
		
		//Only those 12 regions should be kept
		assertTrue(srs.contains(new SimpleRegion("chr1",  21786146 ,  21786776 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  92119247 ,  92119877 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  102969518 ,  102970148 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  121331564 ,  121332194 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  149954581 ,  149955211 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  155787767 ,  155788397 )));
		assertTrue(srs.contains(new SimpleRegion("chr1",  206282567 ,  206283197 )));
		assertTrue(srs.contains(new SimpleRegion("chr12",  100993104 ,  100993734 )));
		assertTrue(srs.contains(new SimpleRegion("chr13",  21729272 ,  21729645 )));
		assertTrue(srs.contains(new SimpleRegion("chr2",  116376661 ,  116377031 )));
		assertTrue(srs.contains(new SimpleRegion("chr4",  169024997 ,  169025627 )));
		
	}
	
	private Vector<MobilePrediction> readInRecords(File fileToRead){
		
		Vector<MobilePrediction> preds = new Vector<MobilePrediction>();
		SAMSilentReader reader = new SAMSilentReader(fileToRead);
		
		for (SAMRecord rec : reader){
			try {
				preds.add(new MobilePrediction(400, 100, 700, rec));
			} catch (IllegalSAMPairException e) {
				e.printStackTrace();
				fail();
			}
		}
		
		reader.close();
		return preds;
	}
	
	
}
