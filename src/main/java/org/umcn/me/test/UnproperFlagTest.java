package org.umcn.me.test;

import net.sf.samtools.*;

import java.io.*;


public class UnproperFlagTest {
	
	public static void main(String[] args){
		validateProperFlag(args[0]);
	
	}
	
	public static void validateProperFlag(String bam){
		SAMFileReader inputSam = new SAMFileReader(new File(bam));
		SAMFileReader inputSam2 = new SAMFileReader(new File(bam));

		int lowest_insert = 1000000;
		int highest_insert = 0;
		
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		inputSam2.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

		
		int c = 0;
		
		for (SAMRecord rec : inputSam){
			c++;
			inputSam2.queryMate(rec);
			System.out.println(c);
//			if (rec.getProperPairFlag()){
//
//				if (!isPairMappedOppositeStrand(rec) ||
//						!isSameReferenceMapped(rec) ||
//						!isFMappedLowerThanR(rec) ||
//						rec.getReadUnmappedFlag() ||
//						rec.getMateUnmappedFlag()){
//					System.out.println("foutje!");
//					System.out.println(rec.getSAMString());
//				}
//				
//				if (Math.abs(rec.getInferredInsertSize()) < lowest_insert ) lowest_insert = rec.getInferredInsertSize();
//				if (Math.abs(rec.getInferredInsertSize()) > highest_insert ) highest_insert = rec.getInferredInsertSize();
//				
//			} else{
//				if (isPairMappedOppositeStrand(rec) && isSameReferenceMapped(rec) && isFMappedLowerThanR(rec)){
//					if (Math.abs(rec.getInferredInsertSize()) >= 116 &&
//							Math.abs(rec.getInferredInsertSize()) <= 70200000){
//						System.out.println(rec.getSAMString());
//					}
//				}
//			}
		}
		System.out.println(lowest_insert);
		System.out.println(highest_insert);
		inputSam.close();

		inputSam2.close();
	}
	
	public static boolean isSameReferenceMapped(SAMRecord read){
		boolean sameReference = false;
		if (read.getReferenceName().equals(read.getMateReferenceName())){
			sameReference = true;
		}
		
		return sameReference;
		
	}
	
	public static boolean isFMappedLowerThanR(SAMRecord read){
		
		boolean mappedLower = false;
		
		if(isPairMappedOppositeStrand(read)){
			if(read.getReadNegativeStrandFlag() && read.getAlignmentStart() >= read.getMateAlignmentStart()){
				mappedLower = true;
			}else if(!read.getReadNegativeStrandFlag() && read.getAlignmentStart() <= read.getMateAlignmentStart()){
				mappedLower = true;
			}
		}

		return mappedLower;
		
	}
	
	public static boolean isPairMappedOppositeStrand(SAMRecord read){
		
		boolean oppStrand = false;
		
		if(!read.getReadNegativeStrandFlag() && read.getMateNegativeStrandFlag()){
			oppStrand = true;
		}else if (read.getReadNegativeStrandFlag() && !read.getMateNegativeStrandFlag()){
			oppStrand = true;
		}
		
		return oppStrand;
	}
}
