package org.umcn.me.tabix;

import java.text.ParseException;

import org.umcn.me.util.SimpleRegion;

public class SelfChainAnnotation extends Annotation {

	public final String qRef;
	public final int qStart;
	public final int qEnd;
	public final double normScore;
	public final String qStrand;
	
	/**
	 * Make a dummy self chain annotation
	 */
	public SelfChainAnnotation(){
		this.qRef = "";
		this.qStart = 0;
		this.qEnd = 0;
		this.normScore = 0.0;
		this.qStrand = "";
	}
	
	public SelfChainAnnotation(String ref, int start, int end, double score, String strand){
		this.qRef = ref;
		this.qStart = start;
		this.qEnd = end;
		this.normScore = score;
		this.qStrand = strand;
	}
	
	@Override
	Annotation parseFromLine(String line) throws ParseException {
		
		String ref;
		int start;
		int end;
		double score;
		
		String strand;
		int size;
		
		String[] split = line.split("\t", -1);
		
		try {
			ref = split[3].trim();
			size = Integer.parseInt(split[4].trim());
			strand = split[5].trim();
			
			int tempStart = Integer.parseInt(split[6].trim());
			int tempEnd = Integer.parseInt(split[7].trim());
			score = Double.parseDouble(split[8].trim());
			
			//TODO: check if calculations are correct
			if (strand.equals("-")){
				start = size - tempEnd + 1; //1-based
				end = size - tempStart; 
			}else{
				start = tempStart + 1; //1-based
				end = tempEnd;
			}
			
		} catch (IndexOutOfBoundsException e) {
			throw new ParseException(e.getMessage(), -1);
		} catch (Exception e){
			throw new ParseException(e.getMessage(), -1);
		}
		
		return new SelfChainAnnotation(ref, start, end , score, strand);
	}
	
	public String toString(){
		return this.qRef + "\t" + this.qStart + "\t" + this.qEnd + "\t" + this.qStrand + "\t" + this.normScore;
	}
	
	public SimpleRegion toRegion(){
		return new SimpleRegion(this.qRef, this.qStart, this.qEnd);
	}
	
	

}
