package org.umcn.me.tabix;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

public class TabixBaseAnnotater<T extends Annotation> {

	protected final TabixReader tr;
	private T dummy;
	
	public TabixBaseAnnotater(TabixReader reader, T dummy){
		this.tr = reader;
		this.dummy = dummy;
	}
	
	public TabixBaseAnnotater(String reader, T dummy) throws IOException{
		this.tr = new TabixReader(reader);
		this.dummy = dummy;
	}
	
	
	@SuppressWarnings("unchecked")
	public List<T> queryOverlapping(String region) throws IOException, ParseException{
		String s;
		List<T> annots = new ArrayList<T>();
		
		TabixReader.Iterator iter = null;
		try {
			iter = tr.query(region);
		} catch (IndexOutOfBoundsException e) {
			System.err.println("Could not query region: " + region);
			return annots;
		}
		
		while (iter != null && (s = iter.next()) != null){
			T annotation = (T) dummy.parseFromLine(s);
			annots.add(annotation);
		}
		
		return annots;
	}

	
}
