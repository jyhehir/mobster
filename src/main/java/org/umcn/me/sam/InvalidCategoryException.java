package org.umcn.me.sam;

public class InvalidCategoryException extends Exception {

	private static final long serialVersionUID = 3124268620353672785L;

	public InvalidCategoryException(String error){
		super(error);
	}

}
