package org.umcn.me.output.vcf;

public class NAConverter {

	public static int toInt(String value){
		if ("NA".equals(value)){
			return 0;
		}
		return Integer.parseInt(value);
	}
	
}
