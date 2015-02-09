package org.umcn.me.util;

import java.io.FileNotFoundException;
import java.security.InvalidParameterException;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import junit.framework.TestCase;

public class BAMCollectionTest extends TestCase {
	
	private final String sam1 = this.getClass().getResource("/test_data/test_sam.sam").getFile();
	private final String sam2 = this.getClass().getResource("/test_data/test_sam2.sam").getFile();
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
		} catch (InvalidParameterException e) {
			fail();
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			fail();
			e.printStackTrace();
		}
	}
	
	//Commented out as the method for now throws an Exception
//	public void testHeaderMergeWhenNoRGPresent(){
//		final String sam1_norg = this.getClass().getResource("/test_data/test_sam_norg.sam").getFile();
//		final String sam2_norg = this.getClass().getResource("/test_data/test_sam2_norg.sam").getFile();
//		final String[] sams_norg = {sam1_norg, sam2_norg};
//		
//		try {
//			BAMCollection col = new BAMCollection(sams_norg, samples);
//			SAMFileHeader header = col.getMergedHeader(SAMFileHeader.SortOrder.unsorted);
//			List<SAMReadGroupRecord> rgs = header.getReadGroups();
//			
//			assertEquals(2, rgs.size());
//			
//		} catch (InvalidParameterException e) {
//			fail();
//			e.printStackTrace();
//		} catch (FileNotFoundException e) {
//			fail();
//			e.printStackTrace();
//		}
//	}
	
}
