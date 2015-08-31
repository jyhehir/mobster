package org.umcn.me.util;

import java.io.File;

public class BAMSample {

	public final File bam;
	public final String sample;
	
	private String prefixReadGroupId = "";
	
	
	public BAMSample (File bam, String sample){
		this.bam = bam;
		this.sample = sample;
	}
	
	/** copy constructor
	 * 
	 * @param bamSample
	 */
	public BAMSample (BAMSample bamSample){
		this.bam = bamSample.bam;
		this.sample = bamSample.sample;
		this.prefixReadGroupId = bamSample.prefixReadGroupId;
	}
	
	public String getPrefixReadGroupId(){
		return this.prefixReadGroupId;
	}
	
	public void setPrefixReadGroupId(final String id){
		this.prefixReadGroupId = id;
	}

	public File getBam(){
		return this.bam;
	}
	
	public String getSample(){
		return this.sample;
	}
	
	public String toString(){
		return this.bam + ":" + this.sample + ":" + this.prefixReadGroupId;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((bam == null) ? 0 : bam.hashCode());
		result = prime * result + ((sample == null) ? 0 : sample.hashCode());
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
		BAMSample other = (BAMSample) obj;
		if (bam == null) {
			if (other.bam != null)
				return false;
		} else if (!bam.equals(other.bam))
			return false;
		if (sample == null) {
			if (other.sample != null)
				return false;
		} else if (!sample.equals(other.sample))
			return false;
		return true;
	}
	
}
