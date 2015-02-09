package org.umcn.me.pairedend;

import java.io.Serializable;

public class TestSerializable implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7718587653359137417L;
	final String name;
	final String city;
	
	public TestSerializable(String s1, String s2){
		this.name = s1;
		this.city = s2;
	}
	
}
