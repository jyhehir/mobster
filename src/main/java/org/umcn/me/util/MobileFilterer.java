package org.umcn.me.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Vector;

import org.apache.log4j.Logger;
import org.umcn.gen.region.*;


public class MobileFilterer {


	public static Logger logger = Logger.getLogger("MobileFilterer");
	
	/**
	 * Create a vector of RegionInterfaces based on the repeatmask file outputted
	 * by UCSC MySQL db
	 * @param repfile path to repeatmask file as outputted by UCSC's MySQL db
	 * @param should this function also sort the contents of the repeatmask file
	 * (use false when the repeatmask file is already sorted)
	 * @return
	 * @throws IOException
	 */
	public static Vector<RegionInterface> createAnnotatedRepeatVector(String repfile, boolean sort) throws IOException{
		
		String line = "";
		String chr = "";
		String[] split;
		int startPos;
		int endPos;

		Vector<RegionInterface> repVector = new Vector<RegionInterface>();
		
		BufferedReader br = new BufferedReader(new FileReader(repfile));
		int c = 0;
		while ((line = br.readLine()) != null){
			c++;
			if(!line.startsWith("bin")){
				split = line.split("\t");
				chr = split[5];
				startPos = Integer.parseInt(split[6].trim()) + 1;
				endPos = Integer.parseInt(split[7].trim());
				
				try {
					repVector.add(new LabeledRegion(Integer.toString(c), chr, startPos, endPos));
				} catch (InvalidRegionException e) {
					logger.warn("Repeat not added to repeatVector, invalid region: " + line);
				}
				

			}
		}

		if(sort){
			Collections.sort(repVector);
		}
		br.close();
		logger.info(repVector.size() + " repeats added to vector");
		return repVector;
	}
		
}
