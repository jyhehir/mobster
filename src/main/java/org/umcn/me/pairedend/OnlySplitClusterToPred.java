package org.umcn.me.pairedend;

import java.io.File;
import java.io.IOException;
import java.util.Vector;

import org.umcn.me.samexternal.IllegalSAMPairException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import static org.umcn.me.output.Save.writePredictionsToFile;

public class OnlySplitClusterToPred {

	public static void main(String[] args){
		SAMFileReader inBAM = new SAMFileReader(new File("C:/bams/BWA/YESSPLIT_r90/ClusterRefacV4_0.1mmYESSplitBWA_splitclustersV2.bam"));
		Vector<MobilePrediction> predictions = new Vector<MobilePrediction>();
		
		for (SAMRecord record : inBAM){
			try {
				MobilePrediction prediction = new MobilePrediction(468, 34, 500, record);
				if (prediction.getLeftTotalHits() + prediction.getRightTotalHits() >= 2){
					predictions.add(prediction);
				}
			} catch (IllegalSAMPairException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		inBAM.close();
		System.out.println("added: " + predictions.size());
		
		try {
			predictions = AnchorClusterer.filterKnownMEs(AnchorClusterer.getKnownMEs(), predictions);
			System.out.println("after filtering: " + predictions.size());
			writePredictionsToFile("C:/ClusterRefacV4_0.1mm5cov_2splithits.txt", predictions, "#testcomment\n", "testsample");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
