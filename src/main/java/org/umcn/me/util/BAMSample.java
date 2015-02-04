package org.umcn.me.util;

import java.io.File;

public class BAMSample {

	public final File bam;
	public final String sample;
	
	public BAMSample (File bam, String sample){
		this.bam = bam;
		this.sample = sample;
	}
	
}
