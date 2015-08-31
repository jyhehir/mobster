package org.umcn.me.parsers;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import org.apache.log4j.Logger;
import org.umcn.gen.sequence.FastaSequence;

/**
 * 
 * Branch of the FastaParserFast class of ExternalInterfaces.
 * Does not need FastaSequence anymore.
 * 
 * @author Djie
 *
 */
public class FastaParserFast {

	public static Logger logger = Logger.getLogger("FastaParserFast");
	
	public static void writeFasta2(Vector<FastaSequence> seqs, String out) throws IOException{
		PrintWriter outFile;
		String header;
		StringBuilder sequence;
		int lineLength = 70;
		int seqLength;
		int end;
		
		outFile = new PrintWriter(new FileWriter(out), true);
		for (FastaSequence seq : seqs){
			header = seq.getHeader().trim();
			sequence = new StringBuilder(seq.getSequence());
			seqLength = sequence.length();
			if (header.startsWith(">")){
				outFile.println(header);
			} else{
				outFile.println(">" + header);
			}
			//Long sequences will be printed as subsequences of lineLength per line
			for(int i = 0; i < seqLength; i += lineLength){
				end = i + lineLength;
				if (end > seqLength){
					end = seqLength;
				}
				outFile.println(sequence.substring(i, end).toUpperCase());
			}			
		}
		outFile.close();
	}
	
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
			FastaParserFast.writeFastQToStream(outFile, rec.getReadName(), sequence, baseQuals);
		}else{
			readId = new StringBuilder(rec.getReadName());
			readId.append(pairDelim);
			if(rec.getFirstOfPairFlag()){
				readId.append("1");
			}else{
				readId.append("2");
			}
			FastaParserFast.writeFastQToStream(outFile, readId.toString(), sequence, baseQuals);
		}
		
	}
}
