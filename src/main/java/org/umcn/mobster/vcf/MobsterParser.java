package org.umcn.mobster.vcf;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.beanio.BeanReader;
import org.beanio.StreamFactory;
import org.beanio.builder.DelimitedParserBuilder;
import org.beanio.builder.StreamBuilder;

public class MobsterParser {

	private File file;
	
	
	public MobsterParser(File mobsterFile){
		this.file = mobsterFile;
	}
	
	
	
	public List<MobsterRecord> parse(){
		
		List<MobsterRecord> records = new ArrayList<MobsterRecord>();
		
		BeanReader in = this.getBeanReader();
		Object record = null;
		in.skip(1); //Skip the header
		while ((record = in.read()) != null){
			MobsterRecord mobsterRec = (MobsterRecord) record;
			records.add(mobsterRec);
		}
		in.close();
		
		return records;
	}
	
	
	private BeanReader getBeanReader(){
		StreamFactory factory = StreamFactory.newInstance();
		
		//Describe the file format
		StreamBuilder builder = new StreamBuilder("mobster")
			.format("delimited")
			.parser(new DelimitedParserBuilder('\t').enableComments("#"))
			.addRecord(MobsterRecord.class); 
		
		factory.define(builder);
		
		BeanReader in = factory.createReader("mobster", this.file);
		
		return in;
	}
	
	
}
