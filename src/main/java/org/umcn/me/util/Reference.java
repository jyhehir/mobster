package org.umcn.me.util;

import net.sf.picard.reference.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

public class Reference {
    private static ReferenceSequenceFile reference = null;
    private static FastaSequenceIndex index = null;

    public static void loadReference(String refPath) throws IOException {

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

        try (BufferedReader refRead = new BufferedReader(new InputStreamReader(new FileInputStream(refPath))); BufferedWriter indexWrite = new BufferedWriter(new FileWriter(indexPath))) {

            //Setup variables
            long offset = 0;
            ArrayList<Integer> lengths = new ArrayList<Integer>();
            int cInt;
            char c;
            boolean append;
            class SequenceIndex {
                String name;
                ArrayList<Integer> lengths;
                long begin;
            }

            //Create indices for each sequence
            SequenceIndex seqIndex = null;
            ArrayList<SequenceIndex> seqIndices = new ArrayList<SequenceIndex>();
            while ((cInt = refRead.read()) != -1) {
                c = (char) cInt;

                //When header...
                if (c == '>') {
                    seqIndex = new SequenceIndex();
                    seqIndex.lengths = new ArrayList<>();
                    offset++;
                    StringBuilder name = new StringBuilder();
                    append = true;

                    //Read it until new line
                    while ((cInt = refRead.read()) != -1) {
                        c = (char) cInt;
                        offset++;
                        if (c == '\n') break;

                        else if (Character.isWhitespace(c)) append = false;
                        if (append) name.append(c);
                    }

                    //Add info to index row
                    seqIndex.name = name.toString();
                    seqIndex.begin = offset;
                    seqIndices.add(seqIndex);

                    System.out.println(name);
                }
                //When sequence...
                else {

                    ///Read it until new line/end of file
                    Integer length = 0;
                    do{
                        c = (char) cInt;
                        offset++;
                        if (c == '\n') break;
                        length++;
                    }
                    while ((cInt = refRead.read()) != -1);

                    //And add info to index row
                    seqIndex.lengths.add(length);
                }
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

    public static String getBaseAt(String chr, long position) throws InvalidNucleotideSequenceException {
        return getSubSequenceAt(chr, position, position).getSequence();
    }

    public static NucleotideSequence getSubSequenceAt(String chr, long start, long stop) throws InvalidNucleotideSequenceException {
        return new NucleotideSequence(new String(reference.getSubsequenceAt(chr, start, stop).getBases(), StandardCharsets.UTF_8));
    }
}
