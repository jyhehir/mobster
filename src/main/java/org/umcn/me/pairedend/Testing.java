package org.umcn.me.pairedend;

import java.io.File;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import net.sf.samtools.*;

import org.umcn.me.util.CollectionUtil;

public class Testing {
	public static void main(String[] args) {
		
		Map<String, String> map = new HashMap<String, String>();
		
		map.put("a", "a");
		map.put("b", "b");
		
		messWithMap(map);
		
		System.out.println(map);
		
		int a = 5;
		int b = -5;
		
		int a_r = - Math.abs(a);
		int b_r = - Math.abs(b);
		
		System.out.println(a_r);
		System.out.println(b_r);
	}

	private static void messWithMap(Map<String, String> map) {
		//map.put("b", "a");
		
		String s1 = map.get("b");

	}
}
