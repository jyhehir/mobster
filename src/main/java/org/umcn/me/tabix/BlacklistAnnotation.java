package org.umcn.me.tabix;

import java.text.ParseException;

public class BlacklistAnnotation extends Annotation {

	public final String component;
	
	/**
	 * Construct dummy annotation
	 */
	public BlacklistAnnotation(){
		this.component = "";
	}
	
	public BlacklistAnnotation(String comp){
		this.component = comp;
	}
	
	@Override
	Annotation parseFromLine(String line) throws ParseException {
		String[] split = line.split("\t", -1);
		String comp;
		try {
			comp = split[3];
		} catch (IndexOutOfBoundsException e) {
			throw new ParseException(e.getMessage(), -1);
		}
		return new BlacklistAnnotation(comp);
	}
	
	public String toString(){
		return this.component;
	}

}
