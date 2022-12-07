package org.umcn.me.output;

import org.apache.log4j.Logger;
import org.umcn.me.output.vcf.MobsterParser;
import org.umcn.me.output.vcf.MobsterRecord;
import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.util.FileValidation;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.SimpleRegion;

import java.io.*;
import java.util.*;

public class FilterPredictions {

    public static Logger logger = Logger.getLogger(FilterPredictions.class.getName());
    private boolean multiple_sample_calling = false;
    private boolean filter_by_read_counts_single_sample = false;
    private static final int KNOWN_ME_FILTER_REGION = 90;
    private static String repmask_file = "../resources/hg19_alul1svaerv.txt";
    private Vector<MobilePrediction> predictions;
    private Map<String, HashMap<String, HashSet<String>>> knownMEs;
    public Map<String, List<MobsterRecord>> normalMEs;
    private String normalPredictionsFile;

    public FilterPredictions(Vector<MobilePrediction> predictions){
        this.predictions = predictions;
    }

    public FilterPredictions(Properties props, Vector<MobilePrediction> predictions) throws IOException {
        if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING) && "true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING).toLowerCase())){
            multiple_sample_calling = true;
        }
        if (props.containsKey(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT)  &&
                "true".equals(props.getProperty(MobileDefinitions.MULTIPLE_SAMPLE_CALLING_STRINGENT).toLowerCase())){
            filter_by_read_counts_single_sample = true;
        }
        this.predictions = predictions;
    }

    public Vector<MobilePrediction> getPredictions(){
        return this.predictions;
    }

    //----- Normal filter methods --------

    public void filterByMinTotalHits(int min_total_hits) {

        int c = 0;
        boolean pass;
        Vector<MobilePrediction> copy_preds = new Vector<MobilePrediction>();
        for (MobilePrediction pred : predictions){

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
        this.predictions = copy_preds;
    }


    /**
     * This function filters predictions when the prediction window overlaps with an annotated
     * mobile element stored in the parameter knownMEs.
     * @param knownMEs - Chromosomes as key
     *                 - Dictionary as value: with mobile element super family as key
     *                        and position start (0-based inclusive) and position end (0-based exclusive)
     *                        as sting value. Start and end are seperated by -
     * @param predictions Vector containing vector of predictions. Each prediction is a string with on
     * index 0: reference chromosome for MEI
     * index 1: mobile family for MEI
     * index 2: estimated insert point
     * index 3: border5 prediction window
     * index 4: border3 prediction window
     * index 5: cluster length of forward cluster (NA if not present)
     * index 6: cluster length of reverse cluster (NA if not present)
     * index 7: number of supporting reads for forward cluster (NA if not present)
     * index 8: number of supporting reads for reverse cluster (NA if not present)
     * @return prediction vector of the same format as arg predictions. Only now with predictions
     * that do not overlap with elements kept in argument knownMEs.
     */

    public void filterKnownMEs(Properties props) throws IOException {
        loadKnownMEs(props);
        filterKnownMEs();
    }

    public void filterKnownMEs(){
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
            border5 = prediction.getLeftPredictionBorder(KNOWN_ME_FILTER_REGION);
            border3 = prediction.getEndRightPredictionBorder(KNOWN_ME_FILTER_REGION);

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
        this.predictions = filteredPredictions;
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
    public void loadKnownMEs(Properties props) throws IOException{
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

        if (!FileValidation.fileValid(repmask_file))
            logger.equals(repmask_file + ": does not exist, or is not readable");
        File f = new File(repmask_file);
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

        this.knownMEs = meMap;
    }

    public void loadNormalMEs(Properties props) throws IOException {

        this.normalPredictionsFile = props.getProperty(MobileDefinitions.NORMAL_PREDICTIONS);

        MobsterParser parser = new MobsterParser(new File(normalPredictionsFile));
        List<MobsterRecord> records = parser.parse();
        Collections.sort(records);

        this.normalMEs = new HashMap<String, List<MobsterRecord>>();
        for(MobsterRecord record: records){
            String chr = record.getChromosome();
            if(!this.normalMEs.containsKey(chr))
                this.normalMEs.put(chr, new ArrayList<MobsterRecord>());
            this.normalMEs.get(chr).add(record);
        }
    }

    public void filterByNormalMEs(Properties props) throws IOException {
        loadNormalMEs(props);
        filterByNormalMEs();
    }

    public void filterByNormalMEs(){

        String ref;
        //String me;
        int border5;
        int border3;
        int alu = 0;
        int sva = 0;
        int l1 = 0;
        int herv = 0;
        int startNormal;
        int endNormal;
        boolean inNormalMe = false;

        Vector<MobilePrediction> filteredPredictions = new Vector<MobilePrediction>();

        for(MobilePrediction prediction : predictions){
            ref = prediction.getOriginalReference();
            border5 = prediction.getLeftPredictionBorder(KNOWN_ME_FILTER_REGION);
            border3 = prediction.getEndRightPredictionBorder(KNOWN_ME_FILTER_REGION);

            if(border5 > border3){
                logger.warn("Safety border5 is bigger than safety border3!");
            }
            for (String me : prediction.getMobileMappings()) {
                if (normalMEs.containsKey(ref)) {
                    for (MobsterRecord normalME : normalMEs.get(ref)) {
                        startNormal = normalME.getBorder5();
                        endNormal = normalME.getEndBorder3();
                        if ((border5 >= startNormal && border5 <= endNormal)
                                || (border3 >= startNormal && border3 <= endNormal)
                                || (border5 <= startNormal && border3 >= endNormal)) {
                            inNormalMe = true;
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
                    if(inNormalMe){
                        break;
                    }
                }
            }
            if(!inNormalMe){
                filteredPredictions.add(prediction);
            } else{
                inNormalMe = false;
            }
        }

        logger.info("filtered alu: " + alu);
        logger.info("filtered sva: " + sva);
        logger.info("filtered l1: " + l1);
        logger.info("filtered herv (hervk): " + herv);
        logger.info("that are already present in the normal predictions file: " + normalPredictionsFile);

        this.predictions = filteredPredictions;
    }


    //----- GRIPS filter methods --------

    /**
     * Only return predictions which occur a maximum of x times based on the name of the source genes in the vector.
     * @param max an integer: only predictions with < max same source genes are kept
     * @return a copy of predictions containing only the predictions which adhere to the filtering criteria.
     */
    public void reducePredictionsBasedOnSource(int max){

        Map<String, Integer> sourceGeneCounts = new HashMap<String, Integer>();
        Set<String> sourceGenesToLog = new HashSet<String>();
        Vector<MobilePrediction> predictionsToKeep = new Vector<MobilePrediction>();
        sourceGeneCounts = getSourceGeneCounts(predictions);


        for (int i=0; i < predictions.size(); i++){
            MobilePrediction pred = predictions.get(i);

            for (String sourceGene : pred.getMobileMappings()){
                if (sourceGeneCounts.containsKey(sourceGene)){
                    int actualCount = sourceGeneCounts.get(sourceGene);

                    if (actualCount < max){
                        predictionsToKeep.add(pred);
                    }
                    if (sourceGenesToLog.add(sourceGene)){
                        logger.debug("Source gene: " + sourceGene +  " - occurs: " + actualCount);
                    }
                }
            }
        }

        logger.info("Separate predictions are allowed to have a maximum of: " + max + " same source genes");
        logger.info("Number of predictions filtered because of same source gene: " + (predictions.size() - predictionsToKeep.size()));

        this.predictions = predictionsToKeep;
    }

    public void removeOverlappingPredictions(){

        Vector<MobilePrediction> predictionsCopy = new Vector<MobilePrediction>(predictions);
        Vector<MobilePrediction> predictionsToReturn = new Vector<MobilePrediction>();


        for (MobilePrediction pred : predictions){
            boolean foundOverlap = false;

            SimpleRegion currentPred = pred.predictionWindowToRegion();

            for (MobilePrediction predCopy : predictionsCopy){
                if (!pred.equals(predCopy) && currentPred.hasOverlapInBP(predCopy.predictionWindowToRegion()) > 0){
                    logger.info("Prediction:  " + pred.getChromosome() + ":" + pred.getInsertionEstimate() + "-" + pred.getMobileMappings().toString() +
                            "and prediction: " + predCopy.getChromosome() + ":" + predCopy.getInsertionEstimate() + "-" + predCopy.getMobileMappings().toString()
                            + " are overlapping and removed");
                    foundOverlap = true;
                    break;
                }
            }

            if (! foundOverlap){
                predictionsToReturn.add(pred);
            }
        }

        logger.info("Number of predictions removed because they were overlapping: " + (predictions.size() - predictionsToReturn.size()));

        this.predictions = predictionsToReturn;
    }

    private Map<String, Integer> getSourceGeneCounts(Vector<MobilePrediction> predictions){

        Map<String, Integer> counts = new HashMap<String, Integer>();

        for (MobilePrediction pred : predictions){
            for (String sourceGene : pred.getMobileMappings()){

                if (counts.containsKey(sourceGene)){
                    int newCount = counts.get(sourceGene) + 1;
                    counts.put(sourceGene, newCount);
                }else{
                    counts.put(sourceGene, 1);
                }

            }
        }
        return counts;

    }
}
