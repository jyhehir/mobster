package org.umcn.me.pairedend;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import org.umcn.me.sam.InvalidCategoryException;
import org.umcn.me.sam.MobileSAMTag;
import org.umcn.me.samexternal.UnknownParamException;

public class RefAndMETest {
	
	public static void main(String[] args) {
		File inFile = new File("D:/wrongParsed.sam");
		SAMFileReader reader = new SAMFileReader(inFile);
		
		String blaat = "";
		
		System.out.println("L1".startsWith(blaat));
		
		for (SAMRecord rec : reader){
			MobileSAMTag tag = new MobileSAMTag();
			tag.setHomoPolymer("polyA");
			try {
				tag.build(rec, "mosaik");
				tag.parse(tag.toString());
			} catch (InvalidCategoryException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (UnknownParamException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println(rec.getSAMString());
			//rec.setAttribute(MobileDefinitions.SAM_TAG_MOBILE, RefAndMEPairFinder.getReadsMappingToME(inFile, new Vector<String>(), "bwa").get("SRGAIIX-593_0004:6:7:8091:9610-1").toString());
			//System.out.println(rec.getSAMString());
			System.out.println(tag);
		}

		reader.close();
	}

}
