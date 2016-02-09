package org.umcn.me.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


import junit.framework.TestCase;

public class RegionWithLabelTest extends TestCase {

	public void testCompare(){
		
		//Same chromosome, same start, diff end
		RegionWithLabel a = new RegionWithLabel("1", "chr1", 100, 1000);
		RegionWithLabel b = new RegionWithLabel("2", "chr1", 100, 2000);
		//b should be > a
		assertTrue(b.compareTo(a) > 0);
		
		//B starts earlier and thus should be < a
		a = new RegionWithLabel("1", "chr1", 100, 1000);
		b = new RegionWithLabel("2", "chr1", 50, 150);
		assertTrue(b.compareTo(a) < 0);
		
		//Now b is exactly the same as a
		a = new RegionWithLabel("1", "chr1", 100, 1000);
		b = new RegionWithLabel("2", "chr1", 100, 1000);
		assertTrue(b.compareTo(a) == 0);
		
		//chr3 is > chr1
		a = new RegionWithLabel("1", "chr1", 100, 1000);
		b = new RegionWithLabel("2", "chr3", 100, 1000);
		assertTrue(b.compareTo(a) > 0);
		
		//we do an alphabetical compare so chr20 is < chr3
		a = new RegionWithLabel("1", "chr20", 100, 1000);
		b = new RegionWithLabel("2", "chr3", 100, 1000);
		assertTrue(b.compareTo(a) > 0);
	}
	
	public void testOverlap(){
		List<RegionWithLabel> regions = new ArrayList<RegionWithLabel>();
		
		RegionWithLabel a = new RegionWithLabel("1", "chr1", 100, 1000);
		RegionWithLabel b = new RegionWithLabel("2", "chr1", 100, 2000);
		RegionWithLabel c = new RegionWithLabel("3", "chr1", 50, 2000);
		RegionWithLabel d = new RegionWithLabel("4", "chr3", 50, 2000);
		RegionWithLabel e = new RegionWithLabel("5", "chr10", 50, 2000);
		RegionWithLabel f = new RegionWithLabel("6", "chr20", 50, 2000);
		RegionWithLabel g = new RegionWithLabel("7", "chrM", 50, 2000);
		
		regions.add(g);
		regions.add(f);
		regions.add(e);
		regions.add(d);
		regions.add(c);
		regions.add(b);
		regions.add(a);
		
		for (RegionWithLabel r : regions){
			System.out.println(r);
		}
		
		a.overlapsWithSortedRegionList(regions);
		
		for (RegionWithLabel r : regions){
			System.out.println(r);
		}
	}
	
	private List<RegionWithLabel> createTestRegions() {
		List<RegionWithLabel> regions = new ArrayList<RegionWithLabel>();
		
		regions.add(new RegionWithLabel("1", "chr1", 100, 200));
		regions.add(new RegionWithLabel("1", "chr1", 300, 400));
		regions.add(new RegionWithLabel("1", "chr1", 500, 600));
		regions.add(new RegionWithLabel("1", "chr2", 90, 210));
		regions.add(new RegionWithLabel("1", "chr2", 100, 200));
		regions.add(new RegionWithLabel("1", "chr2", 120, 180));
		regions.add(new RegionWithLabel("1", "chr2", 300, 400));
		regions.add(new RegionWithLabel("1", "chr2", 500, 600));
		regions.add(new RegionWithLabel("1", "chr3", 100, 200));
		regions.add(new RegionWithLabel("1", "chr3", 150, 250));
		regions.add(new RegionWithLabel("1", "chr3", 200, 300));
		
		regions.add(new RegionWithLabel("1", "chr4", 100, 200));
		regions.add(new RegionWithLabel("1", "chr4", 200, 500));
		regions.add(new RegionWithLabel("1", "chr4", 300, 400));
		regions.add(new RegionWithLabel("1", "chr4", 350, 600));
		regions.add(new RegionWithLabel("1", "chr4", 700, 800));

		return regions;
	}
	
	public void testRegionLookup(){
		List<RegionWithLabel> regions = createTestRegions();
	

		RegionWithLabel r;
		List<RegionWithLabel> oregions = null;
		
		
		r = new RegionWithLabel("1", "chr1", 100, 200);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(1 ,oregions.size());
		
		r = new RegionWithLabel("1", "chr1", 150, 550);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(3, oregions.size());
		
		r = new RegionWithLabel("1", "chr3", 150, 150);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(2, oregions.size());
		
		Collections.sort(regions);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(2, oregions.size());
		
		r = new RegionWithLabel("1", "chr1", 50, 99);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(0, oregions.size());
		
		r = new RegionWithLabel("1", "chr1", 190, 220);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(1, oregions.size());
		
		r = new RegionWithLabel("1", "chr1", 220, 290);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(0, oregions.size());
		
		r = new RegionWithLabel("1", "chr3", 700, 900);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(0, oregions.size());
		
		r = new RegionWithLabel("1", "chr2", 140, 160);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(3, oregions.size());

		r = new RegionWithLabel("1", "chr4", 375, 376);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(3, oregions.size());
		
		r = new RegionWithLabel("1", "chr4", 375, 650);;
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(3, oregions.size());
	
		r = new RegionWithLabel("1", "chr4", 375, 750);
		oregions = r.overlapsWithSortedRegionList(regions);
		assertEquals(4, oregions.size());
		
	}
	

	
}
