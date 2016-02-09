package org.umcn.me.util;

import java.util.ArrayList;
import java.util.List;

public class RegionWithLabel implements Comparable<RegionWithLabel> {

	public final String label;
	public final String chr;
	public final int start;
	public final int end;
	
	public RegionWithLabel(String label, String chrom, int start, int end){
		this.label = label;
		this.chr = chrom;
		this.start = start;
		this.end = end;
	}

	public String toString(){
		return this.chr + "\t" + this.start + "\t" + this.end + "\t" + this.label;
	}
	
	public boolean overlapsWith(RegionWithLabel other){
		boolean overlap = false;		
		if (this.chr.equals(other.chr)){
			//http://stackoverflow.com/questions/3269434/whats-the-most-efficient-way-to-test-two-integer-ranges-for-overlap
			if(this.start <= other.end && other.start <= this.end){
				overlap = true;
			}
		}		
		return overlap;
	}
	
	/**
	 * 
	 * @param regions List must be sorted
	 * @return Overlapping regions are returned. If regions were not sorted then return is undefined.
	 */
	public List<RegionWithLabel> overlapsWithSortedRegionList(List<RegionWithLabel> regions){
		//RegionWithLabel dummyRegion = new RegionWithLabel("1", this.chr, Integer.MIN_VALUE, Integer.MIN_VALUE);
		List<RegionWithLabel> overlap = new ArrayList<RegionWithLabel>();
		for (RegionWithLabel otherRegion : regions){
			if (this.overlapsWith(otherRegion)){
				overlap.add(otherRegion);
			}
		}
		
		return overlap;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + end;
		result = prime * result + ((label == null) ? 0 : label.hashCode());
		result = prime * result + start;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		RegionWithLabel other = (RegionWithLabel) obj;
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (end != other.end)
			return false;
		if (label == null) {
			if (other.label != null)
				return false;
		} else if (!label.equals(other.label))
			return false;
		if (start != other.start)
			return false;
		return true;
	}

	@Override
	public int compareTo(RegionWithLabel other) {
		
		int diff;
		
		if (other == null){
			return 1;
		}
		
		if(this.chr.equals(other.chr)){
			if (this.start == other.start){
				diff = this.end - other.end;
			}else{
				diff = this.start - other.start;
			}
		}else{
			diff = this.chr.compareTo(other.chr);
		}
		
		return diff;
		
	}
	
	
	
	
	
}
