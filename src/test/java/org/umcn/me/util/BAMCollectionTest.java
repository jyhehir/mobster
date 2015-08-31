package org.umcn.me.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMReadGroupRecord;
import junit.framework.TestCase;

public class BAMCollectionTest extends TestCase {
	
	private final String sam1 = this.getClass().getResource("/test_data/test_sam.sam").getFile();
	private final String sam2 = this.getClass().getResource("/test_data/test_sam2.sam").getFile();
	private final String sam_norg = this.getClass().getResource("/test_data/test_sam_norg.sam").getFile();
	private final String[] sams = {sam1, sam2};
	private final String[] samples = {"sample1", "sample2"};
	
	
	public void testHeaderMerge(){
		try {
			BAMCollection col = new BAMCollection(sams, samples);
			SAMFileHeader header = col.getMergedHeader(SAMFileHeader.SortOrder.unsorted);
			List<SAMReadGroupRecord> rgs = header.getReadGroups();
			
			//SAM test file 1 contains 4 read groups, SAM Test file 2 contains 1 group. In total 5 read groups should be created
			assertEquals(5, rgs.size());
			
			//This read group should be set to sample2
			assertEquals("sample2", header.getReadGroup("1").getSample());
			
			System.out.println("Test Header Merge without Prefix");
			
			for (SAMReadGroupRecord rg : rgs){
				System.out.println(rg);
				System.out.println(rg.getId());
			}
			
		} catch (FileNotFoundException e) {
			fail();
			e.printStackTrace();
		}
	}
	
	public void testHeaderMergeWithPrefix(){
		try {
			BAMCollection col = new BAMCollection(sams, samples, true); //prefix readgroups
			SAMFileHeader header = col.getMergedHeader(SAMFileHeader.SortOrder.unsorted);
			List<SAMReadGroupRecord> rgs = header.getReadGroups();
			
			//SAM test file 1 contains 4 read groups, SAM Test file 2 contains 1 group. In total 5 read groups should be created
			assertEquals(5, rgs.size());
			
			//This read group should be set to sample2 + it should be prefixed with 1
			assertEquals("sample2", header.getReadGroup("1.1").getSample());
			
			assertEquals("sample1", header.getReadGroup("0.79454").getSample());
			
			System.out.println("Test Header Merge with Prefix");
			for (SAMReadGroupRecord rg : rgs){
				System.out.println(rg);
				System.out.println(rg.getId());
			}
			
		}  catch (FileNotFoundException e) {
			fail();
			e.printStackTrace();
		}
	}
	
	public void testHeaderMergeWithBAMWithoutReadGroupAndPrefix(){
		try {
			
			String[] samFiles = {sam_norg, sam1};
			String[] samples = {"sample_norg", "sample_rg"};
			
			BAMCollection col = new BAMCollection(samFiles, samples, true); // prefix readgroups
			
			SAMFileHeader header = col.getMergedHeader(SortOrder.unsorted);
			List<SAMReadGroupRecord> rgs = header.getReadGroups();
			
			assertEquals(5, rgs.size());
			assertEquals("sample_norg", header.getReadGroup("0.").getSample());
			assertEquals("0.", col.getPrefixReadGroupIdFromBam(new BAMSample(new File(sam_norg), "sample_norg")));
			

			System.out.println("Test Header Merge with BAM Prefix and a BAM with no readgroup");
			for (SAMReadGroupRecord rg : rgs){
				System.out.println(rg);
				System.out.println(rg.getId());
			}
			
			
		} catch (FileNotFoundException e) {
			fail();
			e.printStackTrace();
		}
	}
	
	
}
