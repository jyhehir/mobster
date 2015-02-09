package org.umcn.me.samexternal;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

public class SAMWriting {

	public static final int MAX_RECORDS_IN_RAM = 10000000;
	
	public static void writeSortedSAMorBAM(File in, File out, File tmp_dir, int max_record, SortOrder order) throws FileNotFoundException{
	
		
		if (!tmp_dir.exists() || !tmp_dir.isDirectory()){
			throw new FileNotFoundException("Temporary directory ( " + tmp_dir.toString() + ") does not exist" );
		}
		
		SAMFileReader reader = new SAMSilentReader(in);
		SAMFileHeader header = reader.getFileHeader();
		
		
		
		//header.setSortOrder(order);
		
		//SAMFileWriterFactory factory = new SAMFileWriterFactory();
		//factory.setMaxRecordsInRam(max_record);
		//factory.setTempDirectory(tmp_dir);
		//factory.setCreateIndex(true);
		SAMFileWriter writer = makeSAMWriter(out, header, tmp_dir, max_record, order);
		
		for (SAMRecord rec: reader){
			writer.addAlignment(rec);
		}
		
		reader.close();
		writer.close();
	}

	public static void writeNameSortedSAMorBAM(File in, File out, File tmp_dir, int max_records) throws FileNotFoundException{
		
		writeSortedSAMorBAM(in, out, tmp_dir, max_records, SortOrder.queryname);
	}

	public static void writeNameSortedSAMorBAM(File in, File out, File tmp_dir) throws FileNotFoundException{
		
		writeNameSortedSAMorBAM(in, out, tmp_dir, MAX_RECORDS_IN_RAM);
	}

	public static void writeNameSortedSAMorBAM(File in, File out, int max_records) throws FileNotFoundException{
		
		writeNameSortedSAMorBAM(in, out, new File(System.getProperty("java.io.tmpdir")), max_records);
	}

	public static void writeNameSortedSAMorBAM(File in, File out) throws FileNotFoundException{
				
		writeNameSortedSAMorBAM(in, out, new File(System.getProperty("java.io.tmpdir")), MAX_RECORDS_IN_RAM);
	}
	
	public static void writeCoordinateSortedBAM(File in, File out, File tmp_dir, int max_records) throws FileNotFoundException{
		writeSortedSAMorBAM(in, out, tmp_dir, max_records, SortOrder.coordinate);
	}

	public static void writeCoordinateSortedBAM(File in, File out, File tmp_dir) throws FileNotFoundException{
		writeCoordinateSortedBAM(in, out, tmp_dir, MAX_RECORDS_IN_RAM);
	}

	public static void writeCoordinateSortedBAM(File in, File out, int max_records) throws FileNotFoundException{
		writeCoordinateSortedBAM(in, out, new File(System.getProperty("java.io.tmpdir")), max_records);
	}

	public static void writeCoordinateSortedBAM(File in, File out) throws FileNotFoundException{
		writeCoordinateSortedBAM(in, out, new File(System.getProperty("java.io.tmpdir")));
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header){
		return makeSAMWriter(out, header, new File(System.getProperty("java.io.tmpdir")));
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, boolean presorted){
		return makeSAMWriter(out, header, new File(System.getProperty("java.io.tmpdir")), MAX_RECORDS_IN_RAM, null, presorted);
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, int max_records){
		return makeSAMWriter(out, header, new File(System.getProperty("java.io.tmpdir")), max_records);
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, File tmp_dir){
		return makeSAMWriter(out, header, tmp_dir, MAX_RECORDS_IN_RAM);
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, File tmp_dir, int max_records){
		return makeSAMWriter(out, header, tmp_dir, max_records, null, false);
	}
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, File tmp_dir, int max_records,
			SortOrder order){
		return makeSAMWriter(out, header, tmp_dir, max_records, order, false);
	}
	
	public static SAMFileWriter makeSAMWriter(File out, SAMFileHeader header, File tmp_dir, int max_records,
			SortOrder order, boolean presorted){
		boolean index = false;
		
		if (order != null){
			header.setSortOrder(order);
		}
		
		if (header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)){
			index = true;
		}
		
		return new SAMFileWriterFactory().setTempDirectory(tmp_dir).setCreateIndex(index).makeBAMWriter(header,
				presorted, out);
	}
}
