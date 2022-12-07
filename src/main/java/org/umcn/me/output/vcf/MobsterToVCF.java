package org.umcn.me.output.vcf;

import java.io.BufferedWriter;
import java.io.File;


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.converters.FileConverter;
import org.umcn.me.util.ReferenceGenome;
import sun.font.TrueTypeFont;

public class MobsterToVCF {

	private String date = "April 8th 2015";
	
	@Parameter(names = "-file", description = "Mobster predictions file (containing metaheader with # and 1 normal line header)", required = true)
	private String in;
	
	@Parameter(names = "-out", description = "Output vcf file", required = true)
	private String out;
	
	public static void main(String[] args) {
		MobsterToVCF mToVCF = new MobsterToVCF();
		JCommander jc = new JCommander(mToVCF);
		jc.setProgramName("MobsterToVCF");
		
		try{ 
			jc.parse(args);
			mToVCF.run();
		} catch (IOException e) {
			System.err.println("MobsterToVCF: could not read in or write outfile.");
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
		run(in, out);
	}



	public static void run(String inString, String outString) throws IOException{
		List<MobsterRecord> records;
		FileWriter fw = new FileWriter(outString);
		BufferedWriter bw = new BufferedWriter(fw);

		MobsterParser parser = new MobsterParser(new File(inString));
		records = parser.parse();
		Collections.sort(records);

		String[] allSamples = new String[0];
		for (MobsterRecord record : records){
			String[] recordSamples = record.getSample().split(", ");
			if(recordSamples.length > allSamples.length){
				allSamples = recordSamples;
			}
		}

		bw.write(MobsterRecordVCFWrapper.getHeader(allSamples));
		for (MobsterRecord record : records){
			MobsterRecordVCFWrapper vcfRecord = new MobsterRecordVCFWrapper(record, allSamples);
			bw.write(vcfRecord.toString());
			bw.write("\n");
		}

		bw.close();
	}

	public static void run(String inString, String outString, ReferenceGenome referenceGenome) throws IOException{
		List<MobsterRecord> records;
		FileWriter fw = new FileWriter(outString);
		BufferedWriter bw = new BufferedWriter(fw);

		MobsterParser parser = new MobsterParser(new File(inString));
		records = parser.parse();
		Collections.sort(records);

		String[] allSamples = new String[0];
		for (MobsterRecord record : records){
			String[] recordSamples = record.getSample().split(", ");
			if(recordSamples.length > allSamples.length){
				allSamples = recordSamples;
			}
		}

		bw.write(MobsterRecordVCFWrapper.getHeader(allSamples));
		for (MobsterRecord record : records){
			MobsterRecordVCFWrapper vcfRecord = new MobsterRecordVCFWrapper(record, allSamples, referenceGenome);
			bw.write(vcfRecord.toString());
			bw.write("\n");
		}

		bw.close();
	}
}
