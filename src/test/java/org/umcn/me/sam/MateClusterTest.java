package org.umcn.me.sam;

import java.io.File;
import java.util.Iterator;
import java.util.Locale;
import java.util.Vector;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import org.umcn.me.samexternal.SAMSilentReader;
import org.umcn.me.util.SampleBam;

import junit.framework.TestCase;

public class MateClusterTest extends TestCase {

	private File anchorFile = new File(this.getClass().getResource("/test_data/A105a_PPIA.sam").getFile());
	
	public void testGetNumberOfDifferentMateMappings(){
		
		//This cluster has mappings to 3 different chromosomes
		MateCluster<SAMRecord> clusterToTest1 = this.anchorsToClusters().get(0);
		
		//This cluster has mappings to 1 chromosome, but also has an unmapped anchor.
		MateCluster<SAMRecord> clusterToTest2 = this.anchorsToClusters().get(2);
		
		//This cluster has mappings to 3 different chromosomes + a couple of unmapped anchors
		MateCluster<SAMRecord> clusterToTest3 = this.anchorsToClusters().get(3);
		
		assertEquals(3, clusterToTest1.getNumberOfDifferentChromosomeMappingsOfMates(true));
		
		assertEquals(1, clusterToTest2.getNumberOfDifferentChromosomeMappingsOfMates(true));
		assertEquals(2, clusterToTest2.getNumberOfDifferentChromosomeMappingsOfMates(false)); //also count the unmapped category
		assertEquals(3, clusterToTest3.getNumberOfDifferentChromosomeMappingsOfMates(true));
		assertEquals(4, clusterToTest3.getNumberOfDifferentChromosomeMappingsOfMates(false));
		
	}
	
	public void testGettingMaxPercentageToSameChromosome(){
		MateCluster<SAMRecord> clusterToTest = this.anchorsToClusters().get(0);
		MateCluster<SAMRecord> clusterToTest3 = this.anchorsToClusters().get(3);
		
		assertEquals("55.56", String.format(Locale.ENGLISH, "%.2f", clusterToTest.getHighestPercentageOfMateAlignmentsToSameChrosome(true)));
		assertEquals(50.0, clusterToTest3.getHighestPercentageOfMateAlignmentsToSameChrosome(true));
		assertEquals(6.0 / 17 * 100, clusterToTest3.getHighestPercentageOfMateAlignmentsToSameChrosome(false));
		
	}
	
	
	private Vector<MateCluster<SAMRecord>> anchorsToClusters(){
		SAMFileReader input = new SAMSilentReader(this.anchorFile);
		SAMFileHeader header = input.getFileHeader();
		SampleBam sampleCalling = SampleBam.SINGLESAMPLE;
		
		
		Vector<MateCluster<SAMRecord>> mobileClusters = new Vector<MateCluster<SAMRecord>>();
		
		double minPercentSameMateRefMapping = 0.0;
		int maxDiffMateMapping = Integer.MAX_VALUE;
		
		Vector<MateCluster<SAMRecord>> clusters = new Vector<MateCluster<SAMRecord>>();
		
		mobileClusters.add(new MateCluster<SAMRecord>(header, false, true, sampleCalling));
		Iterator<MateCluster<SAMRecord>> iter = mobileClusters.iterator();
		
		for(SAMRecord record : input){
			boolean added = false;
			while (iter.hasNext()){
				MateCluster<SAMRecord> currentCluster = iter.next();
				if(currentCluster.isWithinSearchArea(record, 270)){
					if(currentCluster.add(record)){
						added = true;
						break;
					}
				}else if(currentCluster.size() >= 1){
					if (currentCluster.getHighestPercentageOfMateAlignmentsToSameChrosome(true) >= minPercentSameMateRefMapping
							&& currentCluster.getNumberOfDifferentChromosomeMappingsOfMates(true) <= maxDiffMateMapping){
						clusters.add(currentCluster);
					}else{
					}
					iter.remove();
				}else{
					iter.remove();
				}
			}
			
			if (!added){
				MateCluster<SAMRecord> newCluster = new MateCluster<SAMRecord>(header, false, true, sampleCalling);
				newCluster.add(record);
				mobileClusters.add(newCluster);
			}
			iter = mobileClusters.iterator();
		}
		
		input.close();
		for (MateCluster<SAMRecord> cluster : mobileClusters){
			if (cluster.size() >= 1){
				clusters.add(cluster);
			}
		}
		return clusters;
	}
}
