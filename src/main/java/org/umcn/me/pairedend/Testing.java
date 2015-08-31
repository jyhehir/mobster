//package org.umcn.me.pairedend;
//
//import java.io.File;
//import java.util.HashMap;
//
//import org.mapdb.BTreeMap;
//import org.mapdb.DB;
//import org.mapdb.DBMaker;
//import org.mapdb.HTreeMap;
//import org.umcn.gen.sam.SAMSilentReader;
//
//import net.sf.samtools.SAMRecord;
//
//public class Testing {
//	
//	private static final File inBam = new File("C:\\projects_root\\projects\\LTG-lowcov-CNV\\data\\clean\\DNA13-16820-.bam");
//	
//	private static final HashMap<String, String> map = new HashMap<String, String>();
//	
//	public static void main(String[] args) {
//		System.out.println("And Go!");
//		//testRegularHashMap();	
//		testDBHashMap();
//
//	}
//	
//	public static void testRegularHashMap(){
//		
//		SAMSilentReader reader = new SAMSilentReader(inBam);
//		
//		System.out.println("Testing Regular hash map");
//		
//		int c = 0;
//		
//		for (SAMRecord rec : reader){
//			map.put(rec.getReadName(), rec.getSAMString());
//			c++;
//			if (c % 100000 == 0 ){
//				System.out.println(c);
//			}
//		}
//		reader.close();
//	}
//	
//	public static void testDBHashMap(){
//		
//		SAMSilentReader reader = new SAMSilentReader(inBam);
//		
//		System.out.println("Testing mapdb");
//		
//		DB db = DBMaker.newMemoryDB().transactionDisable().make();
//		
//		
//		BTreeMap<String, TestSerializable> m = db.getTreeMap("test");
//		int c = 0;
//		
//		for (SAMRecord rec : reader){
//			TestSerializable test = new TestSerializable("hoi", "blaat");
//			
//			m.put(rec.getReadName(), test);
//			c++;
//			if (c % 100000 == 0 ){
//				System.out.println(c);
//			}
//		}		
//		
//		reader.close();
//	}
//
//}
