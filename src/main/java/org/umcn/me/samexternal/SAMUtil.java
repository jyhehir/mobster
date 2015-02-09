package org.umcn.me.samexternal;

public class SAMUtil {

	
	/**
	 * Converts strings like /data/results/DNA1234049.bam or DNA1234049.bam
	 * to DNA1234049
	 * @param location file path to bam file. String should end in .bam (case insensitive)
	 * @return
	 */
	public static String extractFileNameFromBAMLocation(String location){
		
		String fileName = "";
		location = location.toLowerCase();
		String[] split = location.split("/");
		
		for (String s : split){
			if (s.endsWith(".bam")){
				fileName = s.replaceAll(".bam$", "");
			}
		}
		
		return fileName;
	}
	
}
