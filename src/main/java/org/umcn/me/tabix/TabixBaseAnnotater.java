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
		TabixReader.Iterator iter = tr.query(region);
		String s;
		List<T> annots = new ArrayList<T>();
		while (iter != null && (s = iter.next()) != null){
			T annotation = (T) dummy.parseFromLine(s);
			annots.add(annotation);
		}
		
		return annots;
	}

	
}
