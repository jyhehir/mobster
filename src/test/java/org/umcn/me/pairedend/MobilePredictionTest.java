package org.umcn.me.pairedend;

import junit.framework.TestCase;
import org.junit.Test;
import org.umcn.me.output.FilterPredictions;
import org.umcn.me.util.MobileDefinitions;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.Vector;
import java.util.logging.Filter;

public class MobilePredictionTest extends TestCase {

    // Insertion sites
    // no TSD:
    int no_tsd_left_split = 1000;
    int no_tsd_left_dis = 995;
    int no_tsd_right_split = 1000;
    int no_tsd_right_dis = 1005;
    // Duplication:
    int dup_left_split = 1010;
    int dup_left_dis = 1005;
    int dup_right_split = 990;
    int dup_right_dis = 995;
    // Deletion:
    int del_left_split = 990;
    int del_left_dis = 985;
    int del_right_split = 1010;
    int del_right_dis = 1015;

    //Cluster lengths
    int left_cluster_len = 300;
    int right_cluster_len = 300;

    //Edge case values
    int left_large_cluster_len = 700;
    int right_large_cluster_len = 700;
    int edge_dup_right_dis = 1000;
    int edge_dup_left_dis = 1000;
    int edge2_dup_right_dis = 991;
    int edge2_dup_left_dis = 1009;

    //---Tests for double cluster---

    //Test methods for two split clusters
    @Test
    public void testLsplitRsplit(){
        assertTrue(true);
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                no_tsd_left_split,

                0,
                0,
                0,
                1,
                no_tsd_right_split
        );

        assertEquals(prediction.hasTSD(), "no_tsd");
        assertEquals(no_tsd_left_split, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_left_split, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                dup_left_split,

                0,
                0,
                0,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRsplitDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                del_left_split,

                0,
                0,
                0,
                1,
                del_right_split
        );

        assertEquals("deletion", prediction.hasTSD());
        assertEquals(del_left_split, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_left_split, prediction.getRightPredictionBorder());

        assertEquals(del_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split, prediction.getEndRightPredictionBorder());
    }

    //Test for two split and two discordant clusters
    @Test
    public void testLsplitdisRsplitdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_cluster_len,
                1,
                1,
                no_tsd_left_split,

                no_tsd_right_dis,
                right_cluster_len,
                1,
                1,
                no_tsd_right_split
        );

        assertEquals("no_tsd", prediction.hasTSD());
        assertEquals(no_tsd_left_split, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_left_split, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDisRsplitDisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                dup_left_dis,
                left_cluster_len,
                1,
                1,
                dup_left_split,

                dup_right_dis,
                right_cluster_len,
                1,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDisRsplitDisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                del_left_dis,
                left_cluster_len,
                1,
                1,
                del_left_split,

                del_right_dis,
                right_cluster_len,
                1,
                1,
                del_right_split
        );

        assertEquals("deletion", prediction.hasTSD());
        assertEquals(del_left_split, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_left_split, prediction.getRightPredictionBorder());

        assertEquals(del_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split, prediction.getEndRightPredictionBorder());
    }

    //Test methods for two split and one left discordant cluster
    @Test
    public void testLsplitdisRsplit(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_cluster_len,
                1,
                1,
                no_tsd_left_split,

                0,
                0,
                0,
                1,
                no_tsd_right_split
        );

        assertEquals("no_tsd", prediction.hasTSD());
        assertEquals(no_tsd_left_split, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_left_split, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDisRsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                dup_left_dis,
                left_cluster_len,
                1,
                1,
                dup_left_split,

                0,
                0,
                0,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDisRsplitDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                del_left_dis,
                left_cluster_len,
                1,
                1,
                del_left_split,
                0,
                0,
                0,
                1,
                del_right_split
        );

        assertEquals("deletion", prediction.hasTSD());
        assertEquals(del_left_split, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_left_split, prediction.getRightPredictionBorder());

        assertEquals(del_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split, prediction.getEndRightPredictionBorder());
    }

    //Test methods for two split and one right discordant cluster
    @Test
    public void testLsplitRsplitdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                no_tsd_left_split,

                no_tsd_right_dis,
                right_cluster_len,
                1,
                1,
                no_tsd_right_split
        );

        assertEquals("no_tsd", prediction.hasTSD());
        assertEquals(no_tsd_left_split, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_left_split, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRsplitDisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                dup_left_split,

                dup_right_dis,
                right_cluster_len,
                1,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRsplitDisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                del_left_split,

                del_right_dis,
                right_cluster_len,
                1,
                1,
                del_right_split
        );

        assertEquals("deletion", prediction.hasTSD());
        assertEquals(del_left_split, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_left_split, prediction.getRightPredictionBorder());

        assertEquals(del_right_split, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split, prediction.getEndRightPredictionBorder());
    }

    //Test methods for a left split and right discordant cluster
    @Test
    public void testLsplitRdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                no_tsd_left_split,

                no_tsd_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_split, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_dis, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_dis, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                dup_left_split,

                dup_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_left_split - dup_right_dis < 15? dup_right_dis - (15 - (dup_left_split - dup_right_dis)) :  dup_right_dis, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_dis, prediction.getInsertionEstimate());
        assertEquals(dup_right_dis, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitRdisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                del_left_split,

                del_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_left_split, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_right_dis, prediction.getRightPredictionBorder());

        assertEquals(del_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_dis, prediction.getEndRightPredictionBorder());
    }

    //Test methods for a left discordant and right split cluster
    @Test
    public void testLdisRsplit(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                1,
                no_tsd_right_split
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_dis, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLdisRsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                dup_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(dup_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_dis - dup_right_split < 15? dup_left_dis + (15 - (dup_left_dis - dup_right_split)) :  dup_left_dis, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLdisRsplitDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                del_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                1,
                del_right_split
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_left_dis, prediction.getLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getInsertionEstimate());
        assertEquals(del_right_split, prediction.getRightPredictionBorder());

        assertEquals(del_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split, prediction.getEndRightPredictionBorder());

    }

    //Test methods for two discordant clusters
    @Test
    public void testLdisRdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                no_tsd_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_dis, prediction.getLeftPredictionBorder());
        assertEquals((no_tsd_left_dis+no_tsd_right_dis)/2, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_dis, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals((no_tsd_left_dis+no_tsd_right_dis)/2, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_dis, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLdisRdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                dup_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                dup_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_left_dis - dup_right_dis < 15? dup_right_dis - (15 - (dup_left_dis - dup_right_dis))/2 :  dup_right_dis, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_dis, prediction.getInsertionEstimate());
        assertEquals(dup_right_dis, prediction.getRightPredictionBorder());

        assertEquals(dup_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_dis - dup_right_dis < 15? dup_left_dis + (15 - (dup_left_dis - dup_right_dis))/2 :  dup_left_dis, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLdisRdisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                del_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                del_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_left_dis, prediction.getLeftPredictionBorder());
        assertEquals((del_left_dis+del_right_dis)/2, prediction.getInsertionEstimate());
        assertEquals(del_right_dis, prediction.getRightPredictionBorder());

        assertEquals(del_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals((del_left_dis+del_right_dis)/2, prediction.getEndInsertionEstimate());
        assertEquals(del_right_dis, prediction.getEndRightPredictionBorder());
    }


    //---Test for single clusters---

    //Test methods for one left split cluster
    @Test
    public void testLsplit(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                no_tsd_left_split,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_left_split + 20, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_left_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_left_split + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                dup_left_split,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(dup_left_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getInsertionEstimate());
        assertEquals(dup_left_split + 20, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLsplitDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                del_left_split,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_left_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getInsertionEstimate());
        assertEquals(del_left_split + 20, prediction.getRightPredictionBorder());

        assertEquals(del_left_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(del_left_split, prediction.getEndInsertionEstimate());
        assertEquals(del_left_split + 20, prediction.getEndRightPredictionBorder());

    }

    //Test methods for one left discordant cluster
    @Test
    public void testLdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? no_tsd_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : no_tsd_left_dis, prediction.getInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : no_tsd_left_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? no_tsd_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : no_tsd_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : no_tsd_left_dis + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testLdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                dup_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals(prediction.hasTSD(), "unknown");
        assertEquals(dup_left_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? dup_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : dup_left_dis, prediction.getInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? dup_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : dup_left_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(dup_left_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? dup_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : dup_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? dup_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : dup_left_dis + 20, prediction.getEndRightPredictionBorder());

    }

    @Test
    public void testLdisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                del_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_left_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? del_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : del_left_dis, prediction.getInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? del_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : del_left_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(del_left_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(left_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? del_left_dis + (prediction.median_fragment_length-left_cluster_len+prediction.sd_fragment_length)/2 : del_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(left_cluster_len < prediction.max_expected_cluster_size - 20? del_left_dis + (prediction.max_expected_cluster_size - left_cluster_len) : del_left_dis + 20, prediction.getEndRightPredictionBorder());

    }

    //Test methods for one right split cluster
    @Test
    public void testRsplit(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                0,
                0,
                0,
                1,
                no_tsd_right_split
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_right_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_split + 20, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_right_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(no_tsd_right_split, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_split + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testRsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                0,
                0,
                0,
                1,
                dup_right_split
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(dup_right_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split + 20, prediction.getRightPredictionBorder());

        assertEquals(dup_right_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_right_split + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testRsplitDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                0,
                0,
                0,
                1,
                del_right_split
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(del_right_split - 20, prediction.getLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getInsertionEstimate());
        assertEquals(del_right_split + 20, prediction.getRightPredictionBorder());

        assertEquals(del_right_split - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(del_right_split, prediction.getEndInsertionEstimate());
        assertEquals(del_right_split + 20, prediction.getEndRightPredictionBorder());
    }

    //Test methods for one right discordant cluster
    @Test
    public void testRdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                no_tsd_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_right_dis - (prediction.max_expected_cluster_size - right_cluster_len) : no_tsd_right_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (no_tsd_right_dis + (no_tsd_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : no_tsd_right_dis, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_right_dis - (prediction.max_expected_cluster_size - right_cluster_len) : no_tsd_right_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (no_tsd_right_dis + (no_tsd_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : no_tsd_right_dis, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_dis + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testRdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                dup_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? dup_right_dis - (prediction.max_expected_cluster_size - left_cluster_len) : dup_right_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (dup_right_dis + (dup_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : dup_right_dis, prediction.getInsertionEstimate());
        assertEquals(dup_right_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? dup_right_dis - (prediction.max_expected_cluster_size - left_cluster_len) : dup_right_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (dup_right_dis + (dup_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : dup_right_dis, prediction.getEndInsertionEstimate());
        assertEquals(dup_right_dis + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testRdisDel(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                del_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? del_right_dis - (prediction.max_expected_cluster_size - left_cluster_len) : del_right_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (del_right_dis + (del_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : del_right_dis, prediction.getInsertionEstimate());
        assertEquals(del_right_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(right_cluster_len < prediction.max_expected_cluster_size - 20? del_right_dis - (prediction.max_expected_cluster_size - left_cluster_len) : del_right_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(right_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (del_right_dis + (del_right_dis - (prediction.median_fragment_length - right_cluster_len) - prediction.sd_fragment_length)) / 2 : del_right_dis, prediction.getEndInsertionEstimate());
        assertEquals(del_right_dis + 20, prediction.getEndRightPredictionBorder());
    }

    //Edge case tests
    @Test
    public void testEdgeLdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                no_tsd_left_dis,
                left_large_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(no_tsd_left_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(left_large_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? no_tsd_left_dis + (prediction.median_fragment_length-left_large_cluster_len+prediction.sd_fragment_length)/2 : no_tsd_left_dis, prediction.getInsertionEstimate());
        assertEquals(left_large_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_left_dis + (prediction.max_expected_cluster_size - left_large_cluster_len) : no_tsd_left_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(no_tsd_left_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(left_large_cluster_len < (prediction.median_fragment_length)+ prediction.sd_fragment_length? no_tsd_left_dis + (prediction.median_fragment_length-left_large_cluster_len+prediction.sd_fragment_length)/2 : no_tsd_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(left_large_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_left_dis + (prediction.max_expected_cluster_size - left_large_cluster_len) : no_tsd_left_dis + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testEdgeRdis(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                0,
                0,

                no_tsd_right_dis,
                right_large_cluster_len,
                1,
                0,
                0
        );

        assertEquals("unknown", prediction.hasTSD());
        assertEquals(right_large_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_right_dis - (prediction.max_expected_cluster_size - right_large_cluster_len) : no_tsd_right_dis - 20, prediction.getLeftPredictionBorder());
        assertEquals(right_large_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (no_tsd_right_dis + (no_tsd_right_dis - (prediction.median_fragment_length - right_large_cluster_len) - prediction.sd_fragment_length)) / 2 : no_tsd_right_dis, prediction.getInsertionEstimate());
        assertEquals(no_tsd_right_dis + 20, prediction.getRightPredictionBorder());

        assertEquals(right_large_cluster_len < prediction.max_expected_cluster_size - 20? no_tsd_right_dis - (prediction.max_expected_cluster_size - right_large_cluster_len) : no_tsd_right_dis - 20, prediction.getEndLeftPredictionBorder());
        assertEquals(right_large_cluster_len < prediction.median_fragment_length + prediction.sd_fragment_length? (no_tsd_right_dis + (no_tsd_right_dis - (prediction.median_fragment_length - right_large_cluster_len) - prediction.sd_fragment_length)) / 2 : no_tsd_right_dis, prediction.getEndInsertionEstimate());
        assertEquals(no_tsd_right_dis + 20, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testEdgeLsplitRdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                1,
                dup_left_split,

                edge_dup_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_left_split - edge_dup_right_dis < 15? edge_dup_right_dis - (15 - (dup_left_split - edge_dup_right_dis)) :  edge_dup_right_dis, prediction.getLeftPredictionBorder());
        assertEquals(edge_dup_right_dis, prediction.getInsertionEstimate());
        assertEquals(edge_dup_right_dis, prediction.getRightPredictionBorder());

        assertEquals(dup_left_split, prediction.getEndLeftPredictionBorder());
        assertEquals(dup_left_split, prediction.getEndInsertionEstimate());
        assertEquals(dup_left_split, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testEdgeLdisRsplitDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                edge_dup_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                0,
                0,
                0,
                1,
                dup_right_split
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(dup_right_split, prediction.getLeftPredictionBorder());
        assertEquals(dup_right_split, prediction.getInsertionEstimate());
        assertEquals(dup_right_split, prediction.getRightPredictionBorder());

        assertEquals(edge_dup_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(edge_dup_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(edge_dup_left_dis - dup_right_split < 15? edge_dup_left_dis + (15 - (edge_dup_left_dis - dup_right_split)) :  edge_dup_left_dis, prediction.getEndRightPredictionBorder());
    }

    @Test
    public void testEdgeLdisRdisDup(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                edge2_dup_left_dis,
                left_cluster_len,
                1,
                0,
                0,

                edge2_dup_right_dis,
                right_cluster_len,
                1,
                0,
                0
        );

        assertEquals("duplication", prediction.hasTSD());
        assertEquals(edge2_dup_left_dis - edge2_dup_right_dis < 15? edge2_dup_right_dis - (15 - (edge2_dup_left_dis - edge2_dup_right_dis))/2 :  edge2_dup_right_dis, prediction.getLeftPredictionBorder());
        assertEquals(edge2_dup_right_dis, prediction.getInsertionEstimate());
        assertEquals(edge2_dup_right_dis, prediction.getRightPredictionBorder());

        assertEquals(edge2_dup_left_dis, prediction.getEndLeftPredictionBorder());
        assertEquals(edge2_dup_left_dis, prediction.getEndInsertionEstimate());
        assertEquals(edge2_dup_left_dis - edge2_dup_right_dis < 15? edge2_dup_left_dis + (15 - (edge2_dup_left_dis - edge2_dup_right_dis))/2 :  edge2_dup_left_dis, prediction.getEndRightPredictionBorder());
    }


    @Test
    public void testFilterKnownME(){
        DummyMobilePrediction prediction = new DummyMobilePrediction(
                0,
                0,
                0,
                42,
                153222343,

                0,
                0,
                0,
                0,
                153222307,
                "ALU",
                "chr2"
        );
        Vector<MobilePrediction> predictions = new Vector<>();
        predictions.add(prediction);

        Properties props = new Properties();
        props.put(MobileDefinitions.REPEAT_MASK_FILE, "repmask/alu_l1_herv_sva_other_grch38_ucsc.rpmsk");
        props.put(MobileDefinitions.RESOURCES_DIR, "/Users/ramonvanamerongen/surfdrive/Major research project/programs/mobster2/mobster/resources/");

        FilterPredictions filterPredictions = new FilterPredictions(predictions);
        try {
            filterPredictions.loadKnownMEs(props);
        } catch (IOException e) {
            e.printStackTrace();
        }
        filterPredictions.filterKnownMEs();
        System.out.println(filterPredictions.getPredictions().size());
        assertEquals(0, filterPredictions.getPredictions().size());
    }
}



class DummyMobilePrediction extends MobilePrediction {

    public final static int EXPECTED_DUPLICATION_LENGTH = 15;
    private static Properties props;

    static {
        props = new Properties();
        props.put(MobileDefinitions.SAMPLE_NAME, "test");
    }

    public DummyMobilePrediction(
            int left_mate_cluster_border,
            int left_cluster_length,
            int left_mate_hits,
            int left_aligned_split_hits,
            int left_aligned_split_border,
            int right_mate_cluster_border,
            int right_cluster_length,
            int right_mate_hits,
            int right_aligned_split_hits,
            int right_aligned_split_border,
            String me,
            String chrom
    ) {
        this(left_mate_cluster_border, left_cluster_length, left_mate_hits, left_aligned_split_hits, left_aligned_split_border, right_mate_cluster_border, right_cluster_length, right_mate_hits, right_aligned_split_hits, right_aligned_split_border);
        this.mobile_mappings.add(me);
        this.original_reference = chrom;
    }

    public DummyMobilePrediction(
            int left_mate_cluster_border,
            int left_cluster_length,
            int left_mate_hits,
            int left_aligned_split_hits,
            int left_aligned_split_border,
            int right_mate_cluster_border,
            int right_cluster_length,
            int right_mate_hits,
            int right_aligned_split_hits,
            int right_aligned_split_border
    ) {

        super(props, 470, 35, 500);
        this.left_mate_cluster_border = left_mate_cluster_border;
        this.left_cluster_length = left_cluster_length;
        this.left_mate_hits = left_cluster_length;
        this.left_aligned_split_hits = left_aligned_split_hits;
        this.left_aligned_split_border = left_aligned_split_border;
        this.right_mate_cluster_border = right_mate_cluster_border;
        this.right_cluster_length = right_cluster_length;
        this.right_mate_hits = right_mate_hits;
        this.right_aligned_split_hits = right_aligned_split_hits;
        this.right_aligned_split_border = right_aligned_split_border;
    }
}