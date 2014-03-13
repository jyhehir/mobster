package org.umcn.me.splitread;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class ClippedReadQualDistribution {
	public static void main(String[] args) {
		String[] regions = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
		Vector<String> vectorRegions = new Vector<String>(Arrays.asList(regions));
		File testBam = new File("D:/bams/nipd/Run178_FC2_5500XL_2_20130128_FRAG_BC_DNA1100950.bam");

		SAMFileReader reader = new SAMFileReader(testBam);
		
		Map<Integer, Integer> editBaseStats = new HashMap<Integer, Integer>();
		Map<Integer, Integer> editColorStats = new HashMap<Integer, Integer>();
		Map<Integer, Integer> editHitsStats = new HashMap<Integer, Integer>();
		String editDistanceBase = "NM";
		String colorDistanceBase = "CM";
		String numberHits = "NH";
		
		int max_edit = 0;
		for (SAMRecord record : reader){
			if(record.getCigarString().contains("H")){
				ClippedRead cr;
				try {
					cr = new ClippedRead(record, false, true);
					if (cr.getClippedLength() >= 20){
						System.out.println(cr.getSAMRecord().getSAMString());
					}
				} catch (InvalidHardClipCigarException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(record.getMappingQuality() >= 60){
				int baseEdit = Integer.parseInt(record.getAttribute(editDistanceBase).toString());
				int colorEdit = Integer.parseInt(record.getAttribute(colorDistanceBase).toString());
				int hits = Integer.parseInt(record.getAttribute(numberHits).toString());
//				
//				if(record.getCigarString().contains("H")){
//					ClippedRead cr;
//					try {
//						cr = new ClippedRead(record, false, true);
//						if (cr.getClippedLength() >= 20){
//							System.out.println(cr.getSAMRecord());
//						}
//					} catch (InvalidHardClipCigarException e) {
//						// TODO Auto-generated catch block
//						e.printStackTrace();
//					}
//				}
//				
////				if(editBaseStats.containsKey(baseEdit)){
////					editBaseStats.put(baseEdit, (editBaseStats.get(baseEdit) + 1));
////				}else{
////					editBaseStats.put(baseEdit, 1);
////				}
////				if(editColorStats.containsKey(colorEdit)){
////					editColorStats.put(colorEdit, (editColorStats.get(colorEdit) + 1));
////				}else{
////					editColorStats.put(colorEdit, 1);
////				}
////				if(editHitsStats.containsKey(hits)){
////					editHitsStats.put(hits, (editHitsStats.get(hits) + 1));
////				}else{
////					editHitsStats.put(hits, 1);
////				}
////				
			}
		}
		
//		System.out.println("Base edit distance:");
//		for (int i : editBaseStats.keySet()){
//			System.out.println(i + "\t" + editBaseStats.get(i));
//		}
//		System.out.println("Color edit distance:");
//		for (int i : editColorStats.keySet()){
//			System.out.println(i + "\t" + editColorStats.get(i));
//		}
//		System.out.println("Number hits:");
//		for (int i : editHitsStats.keySet()){
//			System.out.println(i + "\t" + editHitsStats.get(i));
//		}
		//System.out.println(max_edit);
		reader.close();
//		BAMSample testSample = new BAMSample(testBam);
//
//		BAMReadFilter filter = new BAMReadFilter();
//		
//		
//		
//		
//		filter.setMAPQFilter(60);
//		filter.setMaxEditDistance(0);
//		filter.setReferenceRegions(vectorRegions);
//		int totalcounts = 0;
//		Map<String, Integer> counts = testSample.getChromosomeCountsBasedOnFilter(vectorRegions, filter);
//		for (String region : regions){
//			System.out.println(region + ": " + counts.get(region));
//			totalcounts += counts.get(region);
//		}
//		System.out.println("total counts: " + totalcounts);
		
	}
}
