package org.umcn.me.parsers;

import java.io.PrintWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import org.apache.log4j.Logger;

/**
 * 
 * Branch of the FastaParserFast class of ExternalInterfaces.
 * Does not need FastaSequence anymore.
 * 
 * @author Djie
 *
 */
public class FastqParser {

	public static Logger logger = Logger.getLogger("FastaParserFast");
	
	
	//TODO: check whether qual string is a legal qual string
	public static void writeFastQToStream(PrintWriter outFile, String id,
			String seq, String qual){
		
		StringBuilder readId = new StringBuilder();
		
		if (!readId.toString().startsWith("@")){
			readId.append("@");
		}
		
		readId.append(id);
		
		outFile.println(readId.toString());
		outFile.println(seq);
		outFile.println("+");
		outFile.println(qual);
		
	}
	
	public static void writeFastQToStream(PrintWriter outFile, SAMRecord rec, boolean reReverse){
		final String pairDelim = "/"; //default Illumina seperator
		StringBuilder readId;
		String sequence = rec.getReadString();
		String baseQuals = rec.getBaseQualityString();
		
		//make reverse complement of sequence if mapped on negative strand, to get to original fastq record
		if (reReverse && rec.getReadNegativeStrandFlag()){
			sequence = SequenceUtil.reverseComplement(sequence);
			baseQuals = StringUtil.reverseString(baseQuals);
		}
		
		if(!rec.getReadPairedFlag()){
			FastqParser.writeFastQToStream(outFile, rec.getReadName(), sequence, baseQuals);
		}else{
			readId = new StringBuilder(rec.getReadName());
			readId.append(pairDelim);
			if(rec.getFirstOfPairFlag()){
				readId.append("1");
			}else{
				readId.append("2");
			}
			FastqParser.writeFastQToStream(outFile, readId.toString(), sequence, baseQuals);
		}
		
	}
}
