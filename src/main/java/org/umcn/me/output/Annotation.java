package org.umcn.me.output;

import org.umcn.me.tabix.*;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.ReadName;
import org.umcn.me.util.ReadNameOption;
import org.umcn.me.util.ReadNameRetriever;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Properties;

public class Annotation {
    public static List<ReadName> extractReadnames(File anchor, Properties props) {
        boolean usePrefixReference = Boolean.getBoolean(props.getProperty(MobileDefinitions.PREFIX_REFERENCE));
        ReadNameOption option = new ReadNameOption.Builder().addRegion(true).autoPrefixReference(usePrefixReference).build();
        ReadNameRetriever retriever = new ReadNameRetriever(anchor, option);
        List<ReadName> reads = new ArrayList<ReadName>();

        for (ReadName read : retriever){
            reads.add(read);
        }
        return reads;
    }


    //---Annotation methods

    public static void annotateSelfChain(Collection<ReadName> reads, String chainLocation) throws IOException, java.text.ParseException{
        TabixBaseAnnotater<SelfChainAnnotation> tba = new TabixBaseAnnotater<SelfChainAnnotation>(chainLocation, new SelfChainAnnotation());

        for (ReadName name : reads){
            if ( name.mateIsMapped ){
                name.setSelfChain(tba.queryOverlapping(name.toPositionString())); //intentionally query from the anchor position
            }

        }

    }

    public static void annotateRefGene(Collection<ReadName> reads, String transcriptLocation)
            throws IOException, java.text.ParseException {
        TabixBaseAnnotater<RefGeneAnnotation> tba = new TabixBaseAnnotater<RefGeneAnnotation>(transcriptLocation, new RefGeneAnnotation());


        for (ReadName name : reads){
            if( name.mateIsMapped){
                name.setMateRefGeneAnnotation(tba.queryOverlapping(name.mateToPositionString()));
            }
            if (name.isMapped){
                name.setRefGeneAnnotation(tba.queryOverlapping(name.toPositionString()));
            }
        }
    }

    public static void annotateExclusionList(Collection<ReadName> reads, String blacklistLocation) throws IOException, java.text.ParseException{
        TabixBaseAnnotater<BlacklistAnnotation> tba = new TabixBaseAnnotater<BlacklistAnnotation>(blacklistLocation, new BlacklistAnnotation());

        for (ReadName name : reads){
            if (name.isMapped){
                name.setBlacklist(tba.queryOverlapping(name.toPositionString()));
            }
            if (name.mateIsMapped){
                name.setMateBlacklist(tba.queryOverlapping(name.mateToPositionString()));
            }
        }

    }

    public static void annotateRepMask(Collection<ReadName> reads, String repmaskLocation, boolean annotateMates)
            throws IOException, java.text.ParseException {
        TabixBaseAnnotater<RepMaskAnnotation> tba = new TabixBaseAnnotater<RepMaskAnnotation>(repmaskLocation, new RepMaskAnnotation());

        for (ReadName name : reads){
            if(name.isMapped){
                name.setRepMaskAnnotation(tba.queryOverlapping(name.toPositionString()));
            }
            if(annotateMates && name.mateIsMapped){
                name.setMateRepMaskAnnotation(tba.queryOverlapping(name.mateToPositionString()));
            }
        }
    }
}
