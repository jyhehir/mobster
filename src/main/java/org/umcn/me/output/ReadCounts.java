package org.umcn.me.output;

import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.omg.CORBA.DynAnyPackage.InvalidValue;
import org.umcn.me.pairedend.MobilePrediction;
import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.util.ArrayStats;
import org.umcn.me.util.MobileDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class ReadCounts {

    public static Logger logger = Logger.getLogger(ReadCounts.class.getName());
    public static boolean VAF = false;

    public static Vector<MobilePrediction> runFromProperties(Properties props, Vector<MobilePrediction> predictions) throws IOException {

        //Don't predict the VAFs if the index files aren't available
        String samIndexPath = props.getProperty(MobileDefinitions.INFILE_INDEX);
        VAF = samIndexPath != null;

        logger.info("Determining the number of non-supporting reads within the clusters for each prediction");
        if(VAF) logger.info("and calculating the VAF based on the maximum ratio of clipped + discordant / total reads within " + props.getProperty(MobileDefinitions.VAF_DETERMINATION_WINDOW) + " bp to the left and right of the insertion");
        ArrayList<String> samples = new ArrayList<String>(Arrays.asList(props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0)));
        String[] bams = props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0);
        String[] indexes = props.getProperty(MobileDefinitions.INFILE_INDEX).split(MobileDefinitions.DEFAULT_SEP, 0);

        for (MobilePrediction prediction : predictions) {

            //Setup the borders for looking for the VAF of the prediction
            int vafWindow = Integer.parseInt(props.getProperty(MobileDefinitions.VAF_DETERMINATION_WINDOW).trim());
            int leftInsBorder = prediction.getInsertionEstimate();
            int rightInsBorder = prediction.getEndInsertionEstimate() + 1;
            int leftOuterBorder = leftInsBorder - (vafWindow-1);
            int rightOuterBorder = (rightInsBorder+1) + (vafWindow-1);

            //For every sample supporting the prediction, determine the number of non-supporting reads and predict the VAF if possible
            Set<String> predictionSamples = prediction.getSampleNames();
            for(String predictionSample: predictionSamples){
                int samplePos = samples.indexOf(predictionSample);
                String sampleName = samples.get(samplePos);
                String sampleBam = bams[samplePos];
                String sampleIndex = indexes[samplePos];

                //Load the bam file
                File sam = new File(sampleBam);
                File index = new File(sampleIndex);
                SAMSilentReader samFile = new SAMSilentReader(sam, index);

                //Determine the non-supporting reads
                int nonSupportingCount = calculateTotalReadCount(prediction, samFile) - prediction.getSampleCounts().get(sampleName);
                prediction.setNonSupportingCount(sampleName, nonSupportingCount);

                //Calculate the VAF, otherwise set it to '-1'
                if(VAF){
                    double VAF = calculateVAF(prediction.getChromosome(), leftOuterBorder, leftInsBorder,
                            rightInsBorder, rightOuterBorder, samFile);
                    prediction.setVAF(sampleName, VAF);
                } else
                    prediction.setVAF(sampleName, -1);
            }
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

    //Function to calculate the total number of reads that entirely lie in the mate clusters
    //and overlap the clipped positions of the split clusters
    private static int calculateTotalReadCount(MobilePrediction prediction, SAMSilentReader samFile){
        int totalReads = 0;
        if(prediction.getInsertionEstimate() == 125370301){
            System.out.println("R:" + prediction.getRightMateClusterBounderies()[0] + "-" + prediction.getRightMateClusterBounderies()[1]);
        }
        if(prediction.hasLeftMateCluster())
            totalReads += getReadCount(prediction.getChromosome(), prediction.getLeftMateClusterBounderies()[0], prediction.getLeftMateClusterBounderies()[1], samFile, true);
        if(prediction.hasRightMateCluster())
            totalReads += getReadCount(prediction.getChromosome(), prediction.getRightMateClusterBounderies()[0], prediction.getRightMateClusterBounderies()[1], samFile, true);
        if(prediction.hasLeftAlignedSplitCluster() || prediction.hasLeftAlignedSplitCluster())
            totalReads += getReadCount(prediction.getChromosome(), prediction.getSplitClusterBounderies()[0], prediction.getSplitClusterBounderies()[1], samFile, false);
        return totalReads;
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
            ArrayStats.add(totalDepth, readStart, readEnd, 1);

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
                ArrayStats.add(supportDepth, readStart, readEnd, 1);
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
            maxLeftRatio = ArrayStats.max(leftReadRatios);
            if(Double.isNaN(maxLeftRatio) || leftReadRatios.length == 1)
                noLeftRatio = true;
        } catch(InvalidValue e){
            noLeftRatio = true;
        }

        double maxRightRatio = 0;
        boolean noRightRatio = false;

        try{
            maxRightRatio = ArrayStats.max(rightReadRatios);
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

    //Function to retrieve the total number of reads that completely fall within a window
    private static int getReadCount(String referenceSeq, int start, int end, SAMSilentReader samFile, boolean contained){

        //Loop over all reads in the window to count
        int count = 0;
        SAMRecordIterator iter = samFile.query(referenceSeq, start, end, contained);
        while(iter.hasNext()) {
            iter.next();
            count += 1;
        }
        iter.close();

        return count;
    }
}
