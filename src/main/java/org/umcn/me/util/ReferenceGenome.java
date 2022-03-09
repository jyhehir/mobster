package org.umcn.me.util;

import net.sf.picard.reference.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReferenceGenome {
    private static ReferenceSequenceFile reference = null;
    private static FastaSequenceIndex index = null;

    public ReferenceGenome(String refPath) throws IOException {

        //Check if the reference exists
        File refFile = new File(refPath);
        if(!refFile.exists())
            throw new FileNotFoundException("The reference genome was not found");

        //If index doesn't exist at the same location as the reference, create it
        String indexPath = refPath.replaceFirst("\\.[a-z]+$",".fai");

        File indexFile = new File(indexPath);
        if(!indexFile.exists()){
            try {
                createIndex(refPath, indexPath);
            } catch (IOException e) {
                throw new IOException("Could not create reference index file: " + e.toString());
            }
        }

        //Load the index
        FastaSequenceIndex index = new FastaSequenceIndex(indexFile);

        //Load the reference
        reference = new IndexedFastaSequenceFile(new File(refPath), index);
    }

    static public void createIndex(String refPath, String indexPath) throws IOException {

        try (BufferedReader refRead = new BufferedReader(new InputStreamReader(new FileInputStream(refPath)));
             BufferedWriter indexWrite = new BufferedWriter(new FileWriter(indexPath))) {

            //Setup variables
            long offset = 0;
            ArrayList<Integer> lengths = new ArrayList<Integer>();
            String line;
            class SequenceIndex {
                String name;
                ArrayList<Integer> lengths;
                long begin;
            }

            //Create indices for each sequence
            SequenceIndex seqIndex = null;
            ArrayList<SequenceIndex> seqIndices = new ArrayList<SequenceIndex>();
            while ((line = refRead.readLine()) != null) {
                offset += line.length() + 1;

                //When header...
                if (line.startsWith(">")) {
                    seqIndex = new SequenceIndex();
                    seqIndex.lengths = new ArrayList<>();
                    String name;

                    //Extract the name (until the first space)
                    Pattern pattern = Pattern.compile(">([^ ]+)");
                    Matcher matcher = pattern.matcher(line);
                    if(matcher.find())
                        name = matcher.group(1);
                    else
                        throw new IOException("'" + line + "' is not a valid header");

                    //Add info to index row
                    seqIndex.name = name;
                    seqIndex.begin = offset;
                    seqIndices.add(seqIndex);
                }
                //When sequence line, add the length of it to the sequence index
                else
                    seqIndex.lengths.add(line.length());
            }

            //Write the indexes for the sequences to the file
            for(SequenceIndex index: seqIndices) {

                //Check the lengths: They should all be of equal size except for the last one
                if (index.lengths.size() > 1) {
                    for (int i = 1; i < index.lengths.size() - 1; i++) {
                        if (!index.lengths.get(i - 1).equals(index.lengths.get(i)))
                            throw new IOException("Not all sequence lines have equal length");
                    }
                }

                //Write the sequence
                indexWrite.write(index.name + "\t" + index.lengths.stream().reduce(0, Integer::sum) + "\t" + index.begin + "\t" + index.lengths.get(0) + "\t" + (index.lengths.get(0) + 1) + "\n");
            }
        }
    }

    public char getBaseAt(String chr, long position) throws InvalidNucleotideSequenceException {
        return getSubSequenceAt(chr, position, position).getSequence().charAt(0);
    }

    public NucleotideSequence getSubSequenceAt(String chr, long start, long stop) throws InvalidNucleotideSequenceException {
        return new NucleotideSequence(new String(reference.getSubsequenceAt(chr, start, stop).getBases(), StandardCharsets.UTF_8));
    }
}
