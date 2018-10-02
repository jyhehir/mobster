package org.umcn.mobster.vcf;

import java.io.BufferedWriter;
import java.io.File;


import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.converters.FileConverter;

public class MobsterToVCF {

	private String date = "April 8th 2015";
	
	@Parameter(names = "-file", description = "Mobster predictions file (containing metaheader with # and 1 normal line header)", converter = FileConverter.class, required = true)
	private File file;
	
	@Parameter(names = "-out", description = "Output vcf file", converter = FileConverter.class, required = true)
	private File out;
	
	public static void main(String[] args) {
		MobsterToVCF mToVCF = new MobsterToVCF();
		JCommander jc = new JCommander(mToVCF);
		jc.setProgramName("MobsterToVCF");
		
		try{ 
			jc.parse(args);
			mToVCF.run();
		} catch (IOException e) {
			System.err.println("MobsterToVCF: could not write outfile.");
			System.err.println(e.getMessage());
		} catch (ParameterException e){
			if (args.length > 0){
				System.err.println(e.getMessage());
			}
			System.out.println("VERSION: " + mToVCF.date);
			jc.usage();
		}
		
	}

	private void run() throws IOException {
		List<MobsterRecord> records;
		FileWriter fw = new FileWriter(this.out);
		BufferedWriter bw = new BufferedWriter(fw);
		
		MobsterParser parser = new MobsterParser(file);
		records = parser.parse();
		Collections.sort(records);
		
		bw.write(MobsterRecordVCFWrapper.VCFHEADER);
		
		for (MobsterRecord record : records){
			MobsterRecordVCFWrapper vcfRecord = new MobsterRecordVCFWrapper(record);
			bw.write(vcfRecord.toString());
			bw.write("\n");			
		}
		bw.close();
		
	}
	
}
