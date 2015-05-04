package org.umcn.me.util;

/**
 * Simple container class for a reference sequence region.
 * 
 * No checking will be done whether the reference sequence is a known reference sequence.
 * This can be usefull when dealing with references from 
 * 
 * 
 * @author Djie Tjwan Thung
 *
 * TODO: might want to implement RegionInterface... for now not necessary as I use the class only for 
 * region string parsing.
 *
 */
public class SimpleRegion {

	/**
	 * The reference sequence, i.e. chrX
	 */
	public final String reference;
	
	/**
	 * The start position on a reference sequence
	 */
	public final int startPosition;
	
	
	/**
	 * The end position on a reference sequence
	 */
	public final int endPosition;
	
	/**
	 * 
	 * @param region of format <ref>:<start>-<end>
	 */
	public SimpleRegion(String region){
		
		try {
			String[] split = region.split(":");
			String[] split2 = split[1].split("-");
			
			this.reference = split[0];
			this.startPosition = Integer.parseInt(split2[0]);
			this.endPosition = Integer.parseInt(split2[1]);
			
			if (this.endPosition < this.startPosition){
				throw new IllegalArgumentException("The end coordinate should be bigger than the start coordinate");
			}
			
		} catch (Exception e) {
			
			// Catch both NumberFormatException and Index out of bounds
			throw new IllegalArgumentException("Provided region string does not adhere" +
					" to the spec (<ref>:<start>-<end>): " + region);
		}
		
	}
	
	public SimpleRegion(String reference, int start, int end){
		
		if (end < start){
			throw new IllegalArgumentException("The end coordinate should be bigger than the start coordinate");
		}
		if (end < 1 || start < 1){
			throw new IllegalArgumentException("Start and/or end coordinate should not be smaller than 1");
		}
		
		this.reference = reference;
		this.startPosition = start;
		this.endPosition = end;
		
	}
	
	public int getLength(){
		return this.endPosition - this.startPosition + 1;
	}
	
	/**
	 * Print the Simple region in format <ref>:<start>-<end>
	 */
	public String toString(){
		return (this.reference + ":" + this.startPosition + "-" + this.endPosition);
	}
	
	
	public int hasOverlapInBP(SimpleRegion other){
		
		if( ! this.reference.equals(other.reference) ){
			return 0;
		}
		
		/*
		 * Borrowed from RegionComparator
		 * 1: this, 2: other
		 * 
		 * we have a) 1:	---------
		 * 			  2:		----------
		 * 
		 * we have b)  1:	         ---------
		 * 			   2:		----------
		 * 
		 * we have c)  1:	      ------
		 * 			   2:		----------
		 * 
		 * we have d)  1:	    ----------
		 * 			   2:		  ------
		 * 
		 * we have e)  1:	    -----
		 * 			   2:		       ------
		 * 
		 * we have f)  1:	            -----
		 * 			   2:		------
		 * 
		 */
		
		// a
		if (other.startPosition > this.startPosition && other.endPosition > this.endPosition) {
			return (Math.max(0, (this.endPosition - other.startPosition)+1));
		}
		
		// b
		if (this.startPosition > other.startPosition && this.endPosition > other.endPosition) {
			return (Math.max(0, (other.endPosition - this.startPosition)+1));
		}
		
		// c
		if (this.startPosition >= other.startPosition && this.endPosition <= other.endPosition) {
			return this.getLength();
		}
		
		// d
		if (other.startPosition >= this.startPosition && other.endPosition <= this.endPosition) {
			return other.getLength();
		}
		// e, f
		return 0;
		
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + endPosition;
		result = prime * result
				+ ((reference == null) ? 0 : reference.hashCode());
		result = prime * result + startPosition;
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
		SimpleRegion other = (SimpleRegion) obj;
		if (endPosition != other.endPosition)
			return false;
		if (reference == null) {
			if (other.reference != null)
				return false;
		} else if (!reference.equals(other.reference))
			return false;
		if (startPosition != other.startPosition)
			return false;
		return true;
	}
}
