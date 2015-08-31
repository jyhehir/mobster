package org.umcn.me.tabix;

import java.text.ParseException;

public class RepMaskAnnotation extends Annotation {

	public final String repeatName;
	public final String repeatClass;
	public final String repeatFamily;
	
	/**
	 * Construct dummy object
	 */
	public RepMaskAnnotation(){
		this.repeatClass = "";
		this.repeatFamily = "";
		this.repeatName = "";
	}
	
	public RepMaskAnnotation(String repeatName, String repeatClass, String repeatFamily){
		this.repeatName = repeatName;
		this.repeatClass = repeatClass;
		this.repeatFamily = repeatFamily;
	}
	
	@Override
	public RepMaskAnnotation parseFromLine(String line) throws ParseException{
		String[] split = line.split("\t", -1);
		String repName;
		String repClass;
		String repFamily;
		try {
			repName = split[3].trim();
			repClass = split[4].trim();
			repFamily = split[5].trim();
		} catch (IllegalArgumentException e) {
			throw new ParseException(e.getMessage(), -1);
		}
		
		return new RepMaskAnnotation(repName, repClass, repFamily);
	}
	
	public String toString(){
		return this.repeatName + "\t" + this.repeatClass + "\t" + this.repeatFamily;
	}

}
