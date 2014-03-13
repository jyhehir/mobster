package org.umcn.me.pairedend;

import java.io.File;

public class SplitAnchorTest {

	public static void main(String[] args){
		File splitAnchors = new File(args[0]);
		File indexAnchors = new File(args[1]);
		File out = new File(args[2]);

		AnchorClusterer.clusterSplitAnchorsToBAM(splitAnchors, indexAnchors,
				out, 2);
	}
}
