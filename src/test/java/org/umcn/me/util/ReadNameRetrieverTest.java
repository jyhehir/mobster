package org.umcn.me.util;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import org.umcn.me.tabix.TabixReader;
import org.umcn.me.tabix.RefGeneAnnotation;
import org.umcn.me.tabix.RepMaskAnnotation;
import org.umcn.me.tabix.TabixBaseAnnotater;

import junit.framework.TestCase;

public class ReadNameRetrieverTest extends TestCase {

	private final File sam = new File(this.getClass().getResource("/test_data/test_sam.sam").getFile());
	
	public void testRetrievingWhenPrefixingReference(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).prefixReference("chr").build();
		
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		int c = 0;
		for (ReadName name : retriever){
			c++;
			assertTrue(name.reference.startsWith("chr"));
			assertTrue(name.mateReference.startsWith("chr"));
		}
		
		assertEquals(6, c);
	}
	
	public void testRetrievingWhenRemovingPrefixOfReadName(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).prefixLength(2).build();
		
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		int c = 0;
		for (ReadName name : retriever){
			c++;
			if (c == 1){
				assertEquals("IIX-601_0003:5:62:3991:20675", name.readName);
			}
			if (c == 3){
				assertEquals("IIX-601_0003:5:4:14132:5738", name.readName);
			}
		}
		
		assertEquals(6, c);
	}
	
	public void testRetrievingRegions(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).prefixLength(2).build();
		
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		int c = 0;
		
		for (ReadName name : retriever){
			c++;
			if (c == 2){
				assertEquals(229280, name.start);
				assertEquals(229280 + 99, name.end);
				assertEquals(229083, name.mateStart);
				assertEquals(229083 + 99, name.mateEnd);
			}
		}
	}
	
	public void testAutoPrefixingReads(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).autoPrefixReference(true).build();
		
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		int c = 0;
		
		for (ReadName name : retriever){
			c++;
			if (c == 2){
				assertEquals("chr1", name.reference);
			}
			
			if (c == 5){
				assertEquals("chr1", name.reference);
			}
		}
	}
	
	public void testTabixOverlap(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).prefixLength(0).prefixReference("chr").build();
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		try {
			TabixReader tr = new TabixReader("C:/projects_root/projects/Mobster-TCGA/grips/references/hg19_longesttranscripts_27072015_sorted.bed.bgz");
			for (ReadName name : retriever){
				System.out.println(name);
				TabixReader.Iterator iter = tr.query(name.reference + ":" + name.start + "-" + name.end);
				String s;
				while (iter != null && (s = iter.next()) != null){
					System.out.println(s);
				}
				System.out.println("=====");
			}
		} catch (IOException e) {
			fail();
		}

	}
	
	public void testTabixOverlap2(){
		ReadNameOption option = new ReadNameOption.Builder().addRegion(true).prefixLength(0).prefixReference("chr").build();
		ReadNameRetriever retriever = new ReadNameRetriever(this.sam, option);
		try {
			TabixBaseAnnotater<RefGeneAnnotation> tba = new TabixBaseAnnotater<RefGeneAnnotation>("C:/projects_root/projects/Mobster-TCGA/grips/references/hg19_longesttranscripts_27072015_sorted.bed.bgz", new RefGeneAnnotation());
			
			for (ReadName name : retriever){
				List<RefGeneAnnotation> refgenes = tba.queryOverlapping(name.reference + ":" + name.start + "-" + name.end);
				System.out.println(name);
				for (RefGeneAnnotation refgene : refgenes){
					System.out.println(refgene);
				}
				System.out.println("xxxxxx");
//				for (TabixRepMask repmask : repmasks){
//					System.out);
//				}
//				
			}
			
			
		} catch (IOException e) {
			fail();
		} catch (ParseException e) {
			fail();
		}
	}
}
