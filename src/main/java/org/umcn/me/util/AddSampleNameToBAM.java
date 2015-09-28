package org.umcn.me.util;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.umcn.me.samexternal.SAMSilentReader;

public class AddSampleNameToBAM {

	
	public static void main(String[] args) {
		
		HelpFormatter formatter = new HelpFormatter();
		Options options = createOptions();
		String inFile;
		String outFile;
		String sample;
		
		
		if(args.length == 0){
			formatter.printHelp("java -Xmx4g -jar AddSampleNameToBAM.jar" , options);
		}else{
			CommandLineParser parser = new GnuParser();
			
			try {
				CommandLine line = parser.parse(options, args);
				inFile = line.getOptionValue("in");
				sample = line.getOptionValue("sn");
				outFile = line.getOptionValue("out");
				addSampleNameToRecords(inFile, outFile, sample);
				
			} catch (ParseException e) {
				e.printStackTrace();
			}
		}
		
	}
	
	public static void addSampleNameToRecords(String in, String out, String sample){
		SAMFileReader reader = new SAMSilentReader(new File(in));
		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, new File(out));
		
		for (SAMRecord record : reader){
			record.setAttribute(MobileDefinitions.SAM_TAG_SAMPLENAME, sample);
			writer.addAlignment(record);
		}
		
		reader.close();
		writer.close();
	}
	
	private static Options createOptions(){
		Options options = new Options();
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Input BAM file");
		
		options.addOption(OptionBuilder.create("in"));
		
		OptionBuilder.withArgName("String");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("sample name");
		
		options.addOption(OptionBuilder.create("sn"));
		
		OptionBuilder.withArgName("BAM File");
		OptionBuilder.hasArg();
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("Output BAM File");
		
		options.addOption(OptionBuilder.create("out"));
		
		return options;
	}
}
