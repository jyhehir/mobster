package org.umcn.me.tabix;

import java.text.ParseException;

public class RefGeneAnnotation extends Annotation {

	public final String geneSymbol;
	
	/**
	 * Construct dummy object
	 */
	public RefGeneAnnotation(){
		this.geneSymbol = "";
	}
	
	public RefGeneAnnotation(String geneSymbol){
		this.geneSymbol = geneSymbol;
	}
	
	@Override
	Annotation parseFromLine(String line) throws ParseException {
		String[] split = line.split("\t");
		String gene;
		try {
			gene = split[4].trim();
		} catch (IndexOutOfBoundsException e) {
			throw new ParseException(e.getMessage(), -1);
		}
		return new RefGeneAnnotation(gene);
	}
	
	public String toString(){
		return this.geneSymbol;
	}

}
