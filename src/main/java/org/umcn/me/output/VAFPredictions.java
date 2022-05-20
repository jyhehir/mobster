package org.umcn.me.output;

import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.omg.CORBA.DynAnyPackage.InvalidValue;
import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.util.MobileDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Properties;
import java.util.Vector;

public class VAFPredictions {

    public static Logger logger = Logger.getLogger(VAFPredictions.class.getName());

    public static Vector<MobilePrediction> runFromProperties(Properties props, Vector<MobilePrediction> predictions) throws IOException {

        //Don't predict the VAFs when no index file is available
        String samIndexPath = props.getProperty(MobileDefinitions.INFILE_INDEX);
        if(samIndexPath == null)
            return predictions;

        logger.info("Calculating VAF of predictions based on the maximum ratio of clipped + discordant / total reads within " + props.getProperty(MobileDefinitions.VAF_DETERMINATION_WINDOW) + "bp to the left and right of the insertion");
        for (MobilePrediction prediction : predictions) {
            int vafWindow = Integer.parseInt(props.getProperty(MobileDefinitions.VAF_DETERMINATION_WINDOW).trim());

            //Load the bam file
            File sam = new File(props.getProperty(MobileDefinitions.INFILE));
            File index = new File(samIndexPath);
            SAMSilentReader samFile = new SAMSilentReader(sam, index);

            //Setup the borders for looking for the prediction
            int leftInsBorder = prediction.getLeftPredictionBorder();
            int rightInsBorder = prediction.getRightPredictionBorder();

            int leftOuterBorder = leftInsBorder - vafWindow;
            int rightOuterBorder = rightInsBorder + vafWindow;

            //Calculate and setup the VAF
            double VAF = calculateVAF(prediction.getChromosome(), leftOuterBorder, leftInsBorder,
                    rightInsBorder, rightOuterBorder, samFile);
            prediction.setVAF(VAF);
        }

        return predictions;
    }

    //Function to get the VAF based on a region and path to an indexed BAM file
    private static double calculateVAF(String referenceSeq, int leftOuterBorder, int leftInsBorder,
                                       int rightInsBorder, int rightOuterBorder, SAMSilentReader samFile){

        //Determine the ratios of the discordant/split reads and the total depth on the left and right of the insertion point
        double[] leftReadRatios = getReadRatios(referenceSeq, leftOuterBorder, leftInsBorder, samFile);
        double[] rightReadRatios = getReadRatios(referenceSeq, rightInsBorder, rightOuterBorder, samFile);

        //Get the max ratio from both sides
        return getMaxReadRatio(leftReadRatios, rightReadRatios);
    }

    //Function to retrieve the read ratio between the supporting (split/discordant) and total reads
    private static double[] getReadRatios(String referenceSeq, int start, int end, SAMSilentReader samFile){
        int[] totalDepth = new int[end-start+1];
        int[] supportDepth = new int[end-start+1];

        //When the start and the end are the same, return an empty array
        if(end <= start){
            return new double[0];
        }

        //Loop over all reads in the window
        SAMRecordIterator iter = samFile.query(referenceSeq, start, end, false);
        while(iter.hasNext()) {

            //Get the start and end of the read
            SAMRecord read = iter.next();
            int readStart = read.getAlignmentStart();
            int readEnd = read.getAlignmentEnd();

            //Skip any reads that completely fall outside the window
            if( (readStart < start && readEnd < start) ||
                    (readStart > end && readEnd > end)){
                continue;
            }

            //Make the start/end positions of the read relative to the window:
            if(readStart <= start){
                readStart = 0;
            } else {
                readStart = readStart - start;
            }

            if(readEnd >= end){
                readEnd = end-start;
            } else {
                readEnd = readEnd - start;
            }

            //Add the read to the total count
            add(totalDepth, readStart, readEnd, 1);

            //Whenever the read is not a proper pair or it if is clipped add it to the supporting depth
            List<CigarElement> cigar = read.getCigar().getCigarElements();
            boolean isLeftClipped;
            boolean isRightClipped;
            if(cigar.size() > 0){
                isLeftClipped = cigar.get(0).getOperator().equals(CigarOperator.SOFT_CLIP);
                isRightClipped = cigar.get(cigar.size()-1).getOperator().equals(CigarOperator.SOFT_CLIP);
            } else{
                isLeftClipped = false;
                isRightClipped = false;
            }
            if(!read.getProperPairFlag() || isLeftClipped || isRightClipped) {
                add(supportDepth, readStart, readEnd, 1);
            }
        }
        iter.close();

        double[] ratios = new double[end-start+1];
        for(int i=0; i < totalDepth.length; i++){
            ratios[i] = (double)supportDepth[i] / (double)totalDepth[i];
        }

        return ratios;
    }

    //Function to get the maximum ratio of the left and right side
    static double getMaxReadRatio(double[] leftReadRatios, double[] rightReadRatios){
        double maxLeftRatio = 0;
        boolean noLeftRatio = false;

        try{
            maxLeftRatio = max(leftReadRatios);
            if(Double.isNaN(maxLeftRatio) || leftReadRatios.length == 1)
                noLeftRatio = true;
        } catch(InvalidValue e){
            noLeftRatio = true;
        }

        double maxRightRatio = 0;
        boolean noRightRatio = false;

        try{
            maxRightRatio = max(rightReadRatios);
            if(Double.isNaN(maxRightRatio) || rightReadRatios.length == 1)
                noRightRatio = true;
        } catch(InvalidValue e){
            noRightRatio = true;
        }

        //Get the maximum of the two sides as the VAF but set the VAF to -1 when it cannot be determined
        if(noLeftRatio && noRightRatio)
            return -1;
        else if (noLeftRatio)
            return maxRightRatio;
        else if (noRightRatio)
            return maxLeftRatio;
        else
            return Math.max(maxLeftRatio, maxRightRatio);
    }

    //Function to add a number to an array
    static void add(int[] a, int start, int end, int num){
        for(int i = start; i <= end; i++){
            a[i] += num;
        }
    }

    //Function to get the maximum in an array
    static double max(double[] a) throws InvalidValue {

        if(a.length==0){
            throw new InvalidValue("Array is empty");
        }
        double max = a[0];

        for(int i=1; i < a.length; i++){
            if(!Double.isNaN(a[i]) && (a[i] > max || Double.isNaN(max))){
                max = a[i];
            }
        }

        return max;
    }
}
