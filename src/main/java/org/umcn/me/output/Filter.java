package org.umcn.me.output;

import org.apache.log4j.Logger;
import org.umcn.me.pairedend.AnchorClusterer;
import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.util.MobileDefinitions;

import java.io.*;
import java.util.*;

public class Filter {

    public static Logger logger = Logger.getLogger(Filter.class.getName());
    private static boolean multiple_sample_calling = false;
    private static boolean filter_by_read_counts_single_sample = false;
    private static final int FILTER_REGION = 90;
    private static String repmask_file = "./hg19_alul1svaerv.txt";

    public static Vector<MobilePrediction> filterByMinTotalHits(Properties props,
            Vector<MobilePrediction> matePredictions, int min_total_hits) {

        //---Setup variables---
        Vector<MobilePrediction> copy_preds = new Vector<MobilePrediction>();
        int c = 0;
        boolean pass;
        if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING) && "true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING).toLowerCase())){
            multiple_sample_calling = true;
        }
        if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT)  &&
                "true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT).toLowerCase())){
            filter_by_read_counts_single_sample = true;
        }

        for (MobilePrediction pred : matePredictions){

            pass = false;

            if (!multiple_sample_calling || !filter_by_read_counts_single_sample){
                if (pred.getLeftTotalHits() + pred.getRightTotalHits() >= min_total_hits){
                    copy_preds.add(pred);
                    pass = true;
                }
            }

            else{
                Map<String, Integer> sampleCounts = pred.getSampleCounts();
                for (String sample : sampleCounts.keySet()){
                    if (sampleCounts.get(sample) >= min_total_hits){
                        copy_preds.add(pred);
                        pass = true;
                        break;
                    }
                }
            }
            if (!pass){
                c++;
            }
        }

        logger.info(c +" predictions did not meet the minimum required read support of: " + min_total_hits);
        return copy_preds;
    }


    public static Vector<MobilePrediction> filterKnownMEs(Properties props, Vector<MobilePrediction> predictions) throws IOException{
        return filterKnownMEs(getKnownMEs(props), predictions);
    }

    public static Vector<MobilePrediction> filterKnownMEs(Map<String, HashMap<String, HashSet<String>>> knownMEs,
                                                          Vector<MobilePrediction> predictions){
        String ref;
        //String me;
        int border5;
        int border3;
        int alu = 0;
        int sva = 0;
        int l1 = 0;
        int herv = 0;
        int startKnown;
        int endKnown;
        boolean inKnownMe = false;

        Vector<MobilePrediction> filteredPredictions = new Vector<MobilePrediction>();

        for(MobilePrediction prediction : predictions){
            ref = prediction.getOriginalReference();
            //me = prediction.getMateMobileMapping();
            border5 = prediction.getLeftPredictionBorder(FILTER_REGION);
            border3 = prediction.getRightPredictionBorder(FILTER_REGION);

            if(border5 > border3){
                logger.warn("Safety border5 is bigger than safety border3!");
            }

            for (String me : prediction.getMobileMappings()) {
                if (knownMEs.containsKey(ref)
                        && knownMEs.get(ref).containsKey(me)) {
                    for (String coordinates : knownMEs.get(ref).get(me)) {
                        startKnown = Integer
                                .parseInt(coordinates.split("-")[0]);
                        endKnown = Integer.parseInt(coordinates.split("-")[1]);
                        if ((border5 >= startKnown && border5 < endKnown)
                                || (border3 >= startKnown && border3 < endKnown)
                                || (border5 < startKnown && border3 >= endKnown)) {
                            inKnownMe = true;
                            if (me.equals("ALU")) {
                                alu += 1;
                            } else if (me.equals("SVA")) {
                                sva += 1;
                            } else if (me.equals("L1")) {
                                l1 += 1;
                            } else if (me.equals("HERV")) {
                                herv += 1;
                            }
                            break;
                        }
                    }
                    if(inKnownMe){
                        break;
                    }
                }
            }
            if(!inKnownMe){
                filteredPredictions.add(prediction);
            }else{
                inKnownMe = false;
            }
        }
        logger.info("filtered already known alu: " + alu);
        logger.info("filtered already known sva: " + sva);
        logger.info("filtered already known l1: " + l1);
        logger.info("filtered already known herv (hervk): " + herv);
        return filteredPredictions;
    }

    /**
     * This function reads in a RepeatMasker .out file and gets coordinates from the start
     * or end positions (based on index param) of repetitive elements from this file and stores
     * it in a HashMap.
     *
     * @return Map containing chromosomes as keys and HashMap as values, containing
     * mobile elements as keys and position start and position end as string values seperated by "-"
     * @throws IOException
     */
    public static Map<String, HashMap<String, HashSet<String>>> getKnownMEs(Properties props) throws IOException{
        /*
         * Map contains:
         * - Chromosomes as key
         * - Dictionary as value: with mobile element super family as key
         *                        and position start (0-based inclusive) and position end (0-based exclusive)
         *                        as sting value. Start and end are seperated by -
         */
        Map<String, HashMap<String, HashSet<String>>> meMap = new HashMap<String, HashMap<String, HashSet<String>>>();
        InputStream is = null;
        String line = "";
        String chr = "";
        String element = "";
        String[] split;
        String startPos;
        String endPos;
        StringBuilder pos = new StringBuilder();
        boolean header = true;

        if (props.containsKey(MobileDefinitions.REPEAT_MASK_FILE)){
            repmask_file = props.getProperty(MobileDefinitions.RESOURCES_DIR) + props.getProperty(MobileDefinitions.REPEAT_MASK_FILE).trim();
        }

        File f = new File(repmask_file);
        if (!f.exists() || !f.canRead())
            logger.equals(repmask_file + ": does not exist, or is not readable");

        is = new FileInputStream(f);

        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        while ((line = br.readLine()) != null){
            if(!header){
                split = line.split("\t");
                chr = split[5];
                element = split[12].toUpperCase();
                startPos = split[6].trim();
                endPos = split[7].trim();
                pos.append(startPos);
                pos.append("-");
                pos.append(endPos);

                if(element.equals("OTHER")){
                    element = "SVA";
                }else if(element.startsWith("ERV")){
                    element = "HERV";
                }

                if (!meMap.containsKey(chr)){
                    meMap.put(chr, new HashMap<String, HashSet<String>>());
                    meMap.get(chr).put(element, new HashSet<String>());
                    meMap.get(chr).get(element).add(pos.toString());
                }else if(!meMap.get(chr).containsKey(element)){
                    meMap.get(chr).put(element, new HashSet<String>());
                    meMap.get(chr).get(element).add(pos.toString());
                }else{
                    meMap.get(chr).get(element).add(pos.toString());
                }
                pos.setLength(0);
            }
            else{
                header = false;
            }
        }
        br.close();
        is.close();

        return meMap;
    }
}
