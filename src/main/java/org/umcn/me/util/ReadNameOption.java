package org.umcn.me.util;


public class ReadNameOption {

	public final int prefixLength;
	public final boolean addRegion;
	public final String prefixReference;
	
	public ReadNameOption(Builder builder) {
		this.prefixLength = builder.prefixLength;
		this.addRegion = builder.addRegion;
		this.prefixReference = builder.prefixReference;
	}

	public static class Builder{
		
		private int prefixLength = 0;
		private boolean addRegion = true;
		private String prefixReference = "";
		
		public Builder(){
		}
		
		public Builder prefixLength(int length){
			this.prefixLength = length;
			return this;
		}
		
		public Builder addRegion(boolean add){
			this.addRegion = add;
			return this;
		}
		
		public Builder prefixReference(String prefix){
			this.prefixReference = prefix;
			return this;
		}
		
		public ReadNameOption build(){
			return new ReadNameOption(this);
		}
		
	}
	
}
