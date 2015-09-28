package org.umcn.me.util;


public class ReadNameOption {

	public final int prefixLength;
	public final boolean addRegion;
	public final boolean autoPrefixReference;
	public final String prefixReference;
	
	public ReadNameOption(Builder builder) {
		this.prefixLength = builder.prefixLength;
		this.addRegion = builder.addRegion;
		this.prefixReference = builder.prefixReference;
		this.autoPrefixReference = builder.autoPrefixreference;
	}

	public static class Builder{
		
		private int prefixLength = 0;
		private boolean addRegion = true;
		private boolean autoPrefixreference = false;
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
		
		public Builder autoPrefixReference(boolean auto){
			this.autoPrefixreference = auto;
			return this;
		}
		
		public ReadNameOption build(){
			return new ReadNameOption(this);
		}
		
	}
	
}
