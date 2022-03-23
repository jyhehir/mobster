package org.umcn.me.output;

import org.apache.log4j.Logger;
import org.umcn.me.output.vcf.MobsterToVCF;
import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.pairedend.Mobster;
import org.umcn.me.util.CollectionUtil;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.ReadName;

import java.io.*;
import java.util.*;

import static org.umcn.me.output.Annotation.*;

public class SavePredictions {

    public static Logger logger = Logger.getLogger(SavePredictions.class.getName());
    private static final String VERSION;
    private static boolean grips_mode = false;
    private static String outPrefix;
    private static String sample = "";
    private static boolean vcf_out = false;
    private static int min_total_hits = 5;

    static{
        Properties prop = new Properties();
        try {
            prop.load(Mobster.class.getResourceAsStream("/properties/version.properties"));
        } catch (FileNotFoundException e) {
            logger.error("Could not load Mobster properties", e);
        } catch (IOException e) {
            logger.error("Could not read Mobster properties", e);
        }
        VERSION = prop.getProperty("version");
    }

    public static void runFromProperties(Properties props, Vector<MobilePrediction> predictions) throws IOException {

        //---Setup variables---
        if (props.containsKey(MobileDefinitions.GRIPS_MODE) && "true".equals(props.getProperty(MobileDefinitions.GRIPS_MODE).toLowerCase())){
            grips_mode = true;
        }
        if (props.containsKey(MobileDefinitions.MINIMUM_TOTAL_HITS)){
            min_total_hits = Integer.parseInt(props.getProperty(MobileDefinitions.MINIMUM_TOTAL_HITS).trim());
        }
        outPrefix = props.getProperty(MobileDefinitions.OUTFILE).trim();
        if (props.containsKey(MobileDefinitions.VCF) && "true".equals(props.getProperty(MobileDefinitions.VCF).toLowerCase())){
            vcf_out = true;
        }


        //---Start filtering---
        FilterPredictions filter = new FilterPredictions(props, predictions);

        //Remove overlapping predictions
        if (props.containsKey(MobileDefinitions.FILTER_OVERLAPPING_PREDICTIONS) &&
                Boolean.parseBoolean(props.getProperty(MobileDefinitions.FILTER_OVERLAPPING_PREDICTIONS))){
            logger.info("Filter: will remove overlapping predictions");
            filter.removeOverlappingPredictions();
        }

        //When grips is enabled, save unfiltered results
        String commentHeader = getVersionAndParameterInfo(props);
        if (grips_mode){
            logger.info("GRIPS mode enabled, writing unfiltered predictions to: " + outPrefix + "_GRIPS_unfiltered.txt");
            writeGRIPSToFile(outPrefix + "_GRIPS_unfiltered.txt", predictions, commentHeader, false);
        }

        //Remove predictions that have not enough supporting hits
        logger.info("Start filtering...");
        filter.filterByMinTotalHits(min_total_hits);

        //When grips is disabled, remove predictions that overlap known MEIs
        if (!grips_mode){
            filter.filterKnownMEs(props);
        }

        //When grips is enabled, filter for multiple occurring source genes
        if (props.containsKey(MobileDefinitions.GRIPS_MAX_SOURCE_GENES)){
            int maxSourceGenes = Integer.parseInt(props.getProperty(MobileDefinitions.GRIPS_MAX_SOURCE_GENES));
            logger.info("Filtering GRIPS for max source genes: " + maxSourceGenes);
            filter.reducePredictionsBasedOnSource(maxSourceGenes);
        }
        logger.info("End of filtering...");

        predictions = filter.getPredictions();

        //---Save output---
        if (grips_mode){
            logger.info("Starting GRIPS annotation...");

            //Prestep 1 for GRIPS: extracting the read names
            File anchor = new File(props.getProperty(MobileDefinitions.ANCHOR_FILE).trim());
            File splitAnchor = new File(props.getProperty(MobileDefinitions.SPLIT_ANCHOR_FILE).trim());
            List<ReadName> anchorReads = extractReadnames(anchor,props);
            List<ReadName> splitAnchorReads = extractReadnames(splitAnchor,props);
            Map<String, ReadName> readnameMap = new HashMap<String, ReadName>();

            System.out.println("Anchor reads size: " + anchorReads.size());
            System.out.println("Split anchor reads size: " + splitAnchorReads.size());
            readnameMap = CollectionUtil.readNamesToMap(anchorReads);
            System.out.println("read name map size 1: " + readnameMap.size());
            readnameMap.putAll(CollectionUtil.readNamesToMap(splitAnchorReads));
            System.out.println("read name map size 2: " + readnameMap.size());


            //Now only keep the reads in map which made it to the clusters supporting the predictions
            Set<String> reads = new HashSet<String>();
            for (MobilePrediction pred : predictions){
                reads.addAll(pred.getReadNamesFromMateClusters());
                reads.addAll(pred.getReadNamesFromSplitClusters());
            }
            System.out.println("Number of total supporting reads from clusters: " + reads.size());
            Set<String> keys = readnameMap.keySet();
            keys.retainAll(reads);
            System.out.println("Number of keys in map: " + readnameMap.size());

            //annotate the predictions
            String refGeneLoc = props.getProperty(MobileDefinitions.GRIPS_TRANSCRIPT);
            String repMaskLoc = props.getProperty(MobileDefinitions.GRIPS_REPMASK);
            String blackListLoc = props.getProperty(MobileDefinitions.GRIPS_BLACKLIST);
            String selfChainLoc = props.getProperty(MobileDefinitions.GRIPS_SELFCHAIN);

            try {
                annotateRefGene(readnameMap.values(), refGeneLoc);
                annotateRepMask(readnameMap.values(), repMaskLoc, true);
                annotateExclusionList(readnameMap.values(), blackListLoc);
                annotateSelfChain(readnameMap.values(), selfChainLoc);
            } catch (java.text.ParseException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

            //Add in the repmask / refgene annotation
            for (MobilePrediction pred : predictions){
                pred.parseReadNames(readnameMap);
            }
            logger.info("Writing filtered GRIPS predictions to: " + outPrefix + "_GRIPS_predictons.txt and " + outPrefix + "_GRIPS_MORECONFIDENT_predictions.txt");
            writeGRIPSToFile(outPrefix + "_GRIPS_predictons.txt", predictions, commentHeader, false);
            writeGRIPSToFile(outPrefix + "_GRIPS_MORECONFIDENT_predictions.txt", predictions, commentHeader, true);
        }else{
            logger.info("Writing filtered predictions to: " + outPrefix + "_predictions.txt");
            writePredictionsToFile(outPrefix + "_predictions.txt", predictions, commentHeader, sample);
        }

        //When the VCF property is true, convert the predictions file(s) to vcf
        if(vcf_out){
            logger.info("VCF output is enabled, converting generated files to VCF...");
            if(grips_mode) {
                logger.info("Writing unfiltered predictions to: " + outPrefix + "_GRIPS_unfiltered.vcf");
                MobsterToVCF.run(outPrefix + "_GRIPS_unfiltered.txt", outPrefix + "_GRIPS_unfiltered.vcf");
                logger.info("Writing filtered GRIPS predictions to: " + outPrefix + "_GRIPS_predictons.vcf and " + outPrefix + "_GRIPS_MORECONFIDENT_predictions.vcf");
                MobsterToVCF.run(outPrefix + "_GRIPS_predictons.txt", outPrefix + "_GRIPS_predictons.vcf");
                MobsterToVCF.run(outPrefix + "_GRIPS_MORECONFIDENT_predictions.txt", outPrefix + "_GRIPS_MORECONFIDENT_predictions.vcf");
            } else{
                logger.info("Writing filtered predictions to: " + outPrefix + "_predictions.vcf");
                MobsterToVCF.run(outPrefix + "_predictions.txt", outPrefix + "_predictions.vcf");
            }
        }
    }

    public static String getVersionAndParameterInfo(Properties props){
        StringBuilder sb = new StringBuilder();
        Date date = new Date();
        sb.append("#Version: ");
        sb.append(VERSION);
        sb.append("\n");
        sb.append("#Properties file initialization, with following specified values : ");

        for (Object key : props.keySet()){
            String stringKey = (String) key;
            String value = props.getProperty(stringKey);
            sb.append(stringKey);
            sb.append("=");
            sb.append(value);
            sb.append(" ");
        }
        sb.append("\n");
        sb.append("#Creation date: ");
        sb.append(date.toString());
        sb.append("\n");

        return sb.toString();
    }

    public static void writePredictionsToFile(String outString, Vector<MobilePrediction> predictions, String comment, String sampleName){
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(outString), true);
            outFile.print(comment);

            outFile.print(MobilePrediction.getHeader() + "\n");

            for(MobilePrediction pred : predictions){
                //pred.setSampleName(sampleName);
                outFile.print(pred + "\n");
            }

            outFile.close();
        } catch (IOException e) {
            logger.error("Failed to write prediction file: " + e.getMessage());
        }
    }

    public static void writeGRIPSToFile(String outString, Vector<MobilePrediction> predictions, String comment, boolean filter){
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(outString), true);

            boolean writtenHeader = false;
            outFile.print(comment);
            for(MobilePrediction pred : predictions){
                if (! writtenHeader){
                    outFile.println(pred.toGripsHeader());
                    writtenHeader = true;
                }

                if (filter && ! pred.gripNeedsFiltering()){
                    outFile.println(pred.toGripsString());
                }else if (! filter){
                    outFile.println(pred.toGripsString());
                }

            }

            outFile.close();
        } catch (IOException e) {
            logger.error("Failed to write prediction file: " + e.getMessage());
        }
    }

    /**
     * Write a tab delimited file
     *
     * @param lines Nested vector. Each inner vector contains a string value for each column.
     * @param header Vector. Each String inside vector contains a column header.
     * @param outname Name of outputted file.
     * @param delim Delimiter to be used between column names and values.
     * @throws IOException
     */
    public static void writeDelimitedFile(Vector<Vector<String>> lines, Vector<String> header, String outname,
                                          String delim)
            throws IOException{

        PrintWriter outFile = new PrintWriter(new FileWriter(outname), true);
        writeDelimitedLine(header, outFile, delim);

        for(Vector<String> line : lines){
            writeDelimitedLine(line, outFile, delim);
        }

        outFile.close();


    }

    /**
     * Write a single delimited line from a PrintWriter object.
     * @param line Vector containing string values to be printed.
     * @param out PrintWriter object to write from.
     * @param delim Delimiter used for seperating string values in parameter line.
     */
    public static void writeDelimitedLine(Vector<String> line, PrintWriter out, String delim){

        for(int i = 0; i < line.size(); i++){
            out.print(line.get(i));
            if (i != line.size() - 1){
                out.print(delim);
            }
        }
        out.println();
    }
}
