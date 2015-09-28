package org.umcn.me.pairedend;

import java.io.File;
import java.io.IOException;
import java.util.Vector;

import org.umcn.me.samexternal.IllegalSAMPairException;

public class SplitClustersInMateClustersTest {

	public static void main(String[] args){
		
		File clusterIn = new File("D:/bams/test/YESSPLIT_r90/ClusterRefacV4_0.1mmYESSplitBWA_clusters_sorted.bam");
		File clusterIndex = new File("D:/bams/test/YESSPLIT_r90/ClusterRefacV4_0.1mmYESSplitBWA_clusters_sorted.bam.bai");
		File splitCluster = new File("D:/bams/test/YESSPLIT_r90/ClusterRefacV4_0.1mmYESSplitBWA_splitclustersV2sort.bam");
		File splitClusterIn = new File("D:/bams/test/YESSPLIT_r90/ClusterRefacV4_0.1mmYESSplitBWA_splitclustersV2sort.bam.bai");
		
		try {
			Vector<MobilePrediction> predictions = AnchorClusterer.clusterMateClusters(clusterIn, clusterIndex, 50, 700);
			predictions = AnchorClusterer.mergeSplitAndMateClusters(predictions, splitCluster, splitClusterIn);
			predictions = AnchorClusterer.filterKnownMEs(AnchorClusterer.getKnownMEs(), predictions);
			AnchorClusterer.writePredictionsToFile("D:/ClusterRefacV4_0.1mmpSoloSplitSwitchLandRYESSplitBWA_splitmerged.txt", predictions, "#testcomment\n", "testsample");
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalSAMPairException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}

