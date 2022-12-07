package org.umcn.me.pairedend;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.omg.CORBA.DynAnyPackage.InvalidValue;
import org.umcn.me.output.SavePredictions;
import org.umcn.me.output.vcf.MobsterToVCF;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.util.ArrayStats;
import org.umcn.me.util.FileValidation;
import org.umcn.me.util.MobileDefinitions;
import org.umcn.me.util.ReferenceGenome;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.Vector;

public final class Mobster {

    public static Logger logger = Logger.getLogger("Mobster");
    public static String propertiesFile = null;
    public static String inFile = null;
    public static String outFile = null;
    public static String sampleName = null;
    public static String resourcesDir = null;
    public static final String VERSION;

    public static int EXTERNAL_DEPENDENT_PROGRAM_FAIL_EXIT_CODE = -1;
    public static int PROPERTIES_FILE_NOT_FOUND_EXIT_CODE = 1;
    public static int COULD_NOT_READ_PROPERTIES_FILE_EXIT_CODE = 2;

    static {

        //Load version properties
        Properties prop = new Properties();
        try {
            prop.load(Mobster.class.getResourceAsStream("/properties/version.properties"));
        } catch (FileNotFoundException e) {
            logger.error("Could not load version properties", e);
            System.exit(PROPERTIES_FILE_NOT_FOUND_EXIT_CODE);
        } catch (IOException e) {
            logger.error("Could not read version properties", e);
            System.exit(COULD_NOT_READ_PROPERTIES_FILE_EXIT_CODE);
        }
        VERSION = prop.getProperty("version");
    }

    public static void main(String[] args) {

        //Setup configurator
        BasicConfigurator.configure();

        //Setup properties object
        Properties props = new Properties();

        //Parse the arguments and add properties to props
        parseArgs(args, props);

        try {

            //Run picard if enabled
            if (props.containsKey(MobileDefinitions.PICARD_COLLECT_INSERT_METRICS) &&
                    "true".equals(props.getProperty(MobileDefinitions.USE_PICARD).trim())) {
                configureAndRunPicard(props);
            }

            //Extract potential supporting reads
            extractPotentialMEIReads(props);

            //Map them to the mobiome
            final int exitStatus = runMobiomeMapping(props);
            if (exitStatus != 0){
                logger.fatal("Mobiome mapping finished with non-normal exit value: " + exitStatus);
                logger.fatal("Therefore Mobster execution is terminated");
                System.exit(EXTERNAL_DEPENDENT_PROGRAM_FAIL_EXIT_CODE);
            }

            //Extract anchors mapping to mobiome
            runMEIPairFinder(props);

            //Cluster the supporting reads to find the MEI predictions
            Vector<MobilePrediction> predictions = AnchorClusterer.runFromProperties(props);

            //Filter, determine VAF and whether it is somatic/germline (if enabled) and save the predictions
            SavePredictions.runFromProperties(props, predictions);

        } catch (IOException e) {
            logger.error(e.getMessage());
            System.exit(-1);
        }
    }

    // note that this is conditional, i.e. only if user requires this
    private static void configureAndRunPicard(final Properties props) throws IOException {
        String[] sampleBAMs = props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP,0);
        String[] sampleNames = props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP,0);

        //The insert size statistics are averaged over all BAMs
        double[] medians = new double[sampleBAMs.length];
        int[] means = new int[sampleBAMs.length];
        int[] sds = new int[sampleBAMs.length];
        double[] percentiles99 = new double[sampleBAMs.length];

        logger.info("Determining insert size metrics of the sample(s) using picard...");
        for(int i = 0; i < sampleBAMs.length; i++){
            String sampleBAM = sampleBAMs[i].trim();
            String sampleName = sampleNames[i].trim();
            logger.info("Running picard on: " + sampleBAM);
            String picardCommand = "java -Xmx4g -jar " + resourcesDir + props.getProperty(MobileDefinitions.PICARD_COLLECT_INSERT_METRICS) +
                    " VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_" + sampleName + "_hist.pdf" +
                    " INPUT=" + sampleBAM + " OUTPUT=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_" + sampleName + "_insertstats" +
                    " STOP_AFTER=50000000";

            execUnixCommand(picardCommand);

            //Values are on line 8
            //[0 = MEDIAN] [4 = MEAN] [17=99percentile]
            BufferedReader br = new BufferedReader(new FileReader(props.getProperty(MobileDefinitions.OUTFILE).trim() + "_" + sampleName + "_insertstats"));
            String line;
            boolean read = false;
            while ((line = br.readLine()) != null){
                if(read){
                    String[] split = line.split("\t", -1);
                    medians[i] = Double.parseDouble(split[0]);
                    means[i] = (int) Double.parseDouble(split[4]);
                    percentiles99[i] = Double.parseDouble(split[17]);
                    sds[i] = (int) Double.parseDouble(split[5]);
                    break;
                } else if(line.startsWith("MEDIAN_INSERT")){
                    read = true;
                }
            }

            br.close();
        }
        int mean = 0;
        int clustermax = 0;
        int sd = 0;
        double median = 0;
        double percentile99 = 0;
        try{
            mean = ArrayStats.mean(means);
            sd = ArrayStats.mean(sds);
            median = ArrayStats.mean(medians);
            percentile99 = ArrayStats.mean(percentiles99);
            clustermax = (int) (median + percentile99);
        } catch (InvalidValue ignored) {} finally {
            if (mean == 0 || clustermax == 0 || sd == 0){
                logger.error("Could not parse the PICARD CollectInsertSizeMetrics successfully");
                System.exit(-1);
            }
        }

        if (props.containsKey(MobileDefinitions.USE_READ_LENGTH) &&
                "true".equalsIgnoreCase(props.getProperty(MobileDefinitions.USE_READ_LENGTH).trim())){
            clustermax = clustermax - Integer.parseInt(props.getProperty(MobileDefinitions.READ_LENGTH).trim());
        }

        //Only if no valid value was provided, will the picard insert statistics be used
        if(!props.containsKey(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS) || !ArrayStats.isNumeric((String) props.get(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS)))
            props.put(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS, Integer.toString(clustermax));
        if(!props.containsKey(MobileDefinitions.MEAN_FRAGMENT_LENGTH)  || !ArrayStats.isNumeric((String) props.get(MobileDefinitions.MEAN_FRAGMENT_LENGTH)))
            props.put(MobileDefinitions.MEAN_FRAGMENT_LENGTH, Integer.toString(mean));
        if(!props.containsKey(MobileDefinitions.SD_FRAGMENT_LENGTH) || !ArrayStats.isNumeric((String) props.get(MobileDefinitions.SD_FRAGMENT_LENGTH)))
            props.put(MobileDefinitions.SD_FRAGMENT_LENGTH, Integer.toString(sd));

        logger.info("Insert size metrics to be used are MEAN_FRAGMENT_LENGTH=" + props.getProperty(MobileDefinitions.MEAN_FRAGMENT_LENGTH) + ", SD_FRAGMENT_LENGTH=" + props.getProperty(MobileDefinitions.SD_FRAGMENT_LENGTH)  + " and LENGTH_99PROCENT_OF_FRAGMENTS=" + props.getProperty(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS) );
    }

    private static void extractPotentialMEIReads(final Properties props) throws IOException {
        props.put(MobileDefinitions.INFILE_FROM_POTENTIAL_MEI_FINDER,
                props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.bam");
        PotentialMEIReadFinder.runFromProperties(props);
    }

    //do the mobiome mapping here
    private static int runMobiomeMapping(final Properties props) {
        String mobiomeMappingCmd = props.getProperty(MobileDefinitions.MOBIOME_MAPPING_CMD).trim();

        // note that this assumes that upstream task RefAndMEPairFinder outputs "_potential.fq"
        mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(FASTQ\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.fq");

        if (mobiomeMappingCmd.toLowerCase().contains("mosaik")) {
            mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(DAT_FILE\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_potential.dat"); // this line is specific to MOSAIK
            mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(RESOURCES_DIR\\)", resourcesDir == null? "" : resourcesDir); // this line is specific to MOSAIK
        }

        // note that this assumes the aligner used will automatically append ".bam", i.e. it will output "_mappedpotentials.bam"
        mobiomeMappingCmd = mobiomeMappingCmd.replaceAll("\\(OUT_FILE\\)", props.getProperty(MobileDefinitions.OUTFILE).trim() + "_mappedpotentials");

        props.put(MobileDefinitions.INFILE_FROM_MOBIOME_MAPPING,
                props.getProperty(MobileDefinitions.OUTFILE).trim() + "_mappedpotentials.bam");

        return execUnixCommand(mobiomeMappingCmd);
    }

    private static final int CAUGHT_IOEXCEPTION_IN_EXECUTION = -1;
    private static final int PROCESS_INTERRUPTED = -2;
    private static int execUnixCommand(String cmd) {

        ProcessBuilder builder = new ProcessBuilder("/bin/sh", "-c", cmd);
        builder.redirectErrorStream(true);
        try {
            String s;
            Process process = builder.start();
            BufferedReader stdout = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            while ((s = stdout.readLine()) != null) {
                System.out.println(s);
            }
            int exitStatus = process.waitFor();

            System.out.println("Exit value: " + exitStatus);

            return exitStatus;

        } catch (IOException e) {
            logger.error("Error in executing command: " + cmd);
            logger.error(e.getMessage());
            return CAUGHT_IOEXCEPTION_IN_EXECUTION;
        } catch (InterruptedException e) {
            logger.error("Interrupted while executing command: " + cmd);
            logger.error(e.getMessage());
            return PROCESS_INTERRUPTED;
        }
    }

    private static void runMEIPairFinder(final Properties props) throws IOException {
        props.put(MobileDefinitions.ANCHOR_FILE, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_anchors.bam");
        props.put(MobileDefinitions.SPLIT_ANCHOR_FILE, props.getProperty(MobileDefinitions.OUTFILE).trim() + "_splitanchors.bam");
        RefAndMEPairFinder.runFromPropertiesFile(props);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static void printUsage() {
        System.out.println("##########################");
        System.out.println("#MOBSTER                 #");
        System.out.println("##########################");
        System.out.println("Version: " + VERSION);
        System.out.println("Author: Djie Tjwan Thung");
        System.out.println();
        System.out.println("Predict non-reference Mobile Element Insertion (MEI) events using one properties file.");
        System.out.println("Only the properties file is required, but it is advisable to also provide the other single dash ('-[argument]') arguments if applicable.");
        System.out.println("These other arguments will override the default properties in the properties file.");
        System.out.println("In addition, any property as listed in the default properties file can also be overridden by providing it as a double dash argument ('--USE_PICARD [value]').");
        System.out.println("A new properties file containing all the overridden properties will also be created with the output prefix.");
        System.out.println();
        System.out.println("\t-properties [properties]\t\tThis value is required. Path to the properties file.");
        System.out.println(("\t-rdir [resources folder]\t\tSpecify where the resources folder of mobster is installed."));
        System.out.println("\t-in [input .bam file]\t\t\tMultiple BAM files may be specified if separated by a comma.");
        System.out.println("\t-out [output prefix]\t\t\tPrefix for the output files.");
        System.out.println("\t-sn [sample name]\t\t\t\tMultiple sample names may be specified if separated by a comma. When tumor/normal paired samples are given, the normal sample should be specified with '-ns'.");
        System.out.println("\t-r [reference genome]\t\t\tWhen provided, sequences will be determined for the target site duplication and the insertion position.");
        System.out.println("\t-np [normal predictions file]\tCan be used to predict somatic MEIs if only the tumor file is provided instead of providing both a tumor and normal sample. If the file provided (.vcf or .txt ) contains Mobster's predictions of the normal sample, any new prediction in the tumor sample (provided through '-in') near a prediction in this file will be removed.");
        System.out.println("\t-ns [normal sample name]\t\t#Can be used to predict somatic MEIs if both a tumor and normal sample have been provided through '-in' and '-sn'. The sample with the corresponding name among the provided sample names will be seen as the normal sample. The other sample will be seen as the tumor sample.");
        System.out.println("\t-vcf\t\t\t\t\t\t\tIf this flag is given, the results of Mobster are outputted in VCF format in addition to the default format.");
        System.out.println(("\t--<PROPERTY> [value]\t\t\tThis value will override the corresponding value in the properties file."));
        System.out.println();
        System.out.println("Default mapping tool: " + SAMDefinitions.MAPPING_TOOL_UNSPECIFIED);
        System.out.println();
    }

    private static void parseArgs(String[] args, Properties props) {

        try{
            if (args.length > 0 && args[0] != null) {

                //Loop over all properties to find if a properties file has been provided, load it
                //And add the properties to the properties object
                for (int i = 0; i < args.length; i += 1)
                    if (args[i].equalsIgnoreCase("-properties")) {
                        propertiesFile = args[i + 1];
                        props.load(new FileInputStream(propertiesFile));
                    }

                //Again loop over all arguments and override any properties in props when they are provided on the command line
                for (int i = 0; i < args.length; i += 1) {
                    String argument = args[i];

                    //First check for flags
                    if(argument.equalsIgnoreCase("-vcf")){
                        props.put(MobileDefinitions.VCF, "true");
                        continue;
                    }

                    //Then for arguments
                    String argUpperCase = argument.replace("--", "").toUpperCase();
                    String value = args[++i];
                    if (argument.equalsIgnoreCase("-in")){
                        inFile = value;
                        props.put(MobileDefinitions.INFILE, inFile);
                    } else if (argument.equalsIgnoreCase("-out")){
                        outFile = value;
                        props.put(MobileDefinitions.OUTFILE, outFile);
                    } else if (argument.equalsIgnoreCase("-sn")){
                        sampleName = value;
                        props.put(MobileDefinitions.SAMPLE_NAME, sampleName);
                    }
                    else if (argument.equalsIgnoreCase("-r")){
                        props.put(MobileDefinitions.REFERENCE_GENOME_FILE, value);
                    }
                    else if (argument.equalsIgnoreCase("-np")){
                        props.put(MobileDefinitions.NORMAL_PREDICTIONS, value);
                    }
                    else if (argument.equalsIgnoreCase("-rdir")){
                        props.put(MobileDefinitions.RESOURCES_DIR, value);
                    }
                    else if (argument.equalsIgnoreCase("-ns")){
                        props.put(MobileDefinitions.NORMAL_SAMPLE, value);
                    }
                    //When preceded by '--' it is presumed to be a property and it is overridden.
                    else if(argument.startsWith("--")){
                        if(props.containsKey(argUpperCase))
                            logger.info("The property \"" + argUpperCase + "\" has been set to \"" + value + "\" from \"" + props.get(argUpperCase) + "\"");
                        else
                            logger.info("The property \"" + argUpperCase + "\" has been set to \"" + argUpperCase + "\" but was not yet present in the properties file.");
                        props.put(argUpperCase, value);

                    } else if(!argument.equals("-properties")){
                        printUsage();
                        logger.error("Invalid argument: \"" + argument +"\". Please try again.");
                        System.exit(-1);
                    }
                }

                //Write the current properties to a file
                String output = props.getProperty(MobileDefinitions.OUTFILE) + "_Mobster.properties";
                props.store(new FileOutputStream(output), "");

                //Check whether there is a correct number of samples
                if (props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0).length !=
                        props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0).length){
                    printUsage();
                    logger.error("Number of supplied samples does not equal the number of supplied bams. Exiting");
                    System.exit(1);
                } else if (props.getProperty(MobileDefinitions.SAMPLE_NAME).split(MobileDefinitions.DEFAULT_SEP, 0).length > 1){
                    logger.info("Detected multiple samples. Multiple sample calling will be turned ON even if this property was set to false in the properties file");
                    props.put(MobileDefinitions.MULTIPLE_SAMPLE_CALLING, "true");
                }

                //Check if the provided sample BAMs file have an index file
                //If not VAF determination will be disabled
                boolean noBamIndex = false;
                if(props.containsKey(MobileDefinitions.INFILE_INDEX)) {
                    if (props.getProperty(MobileDefinitions.INFILE_INDEX).split(MobileDefinitions.DEFAULT_SEP, 0).length !=
                            props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0).length) {
                        logger.warn("The number of samples does not equal the number of provided BAM index files. VAF determination will be disabled.");
                        props.remove(MobileDefinitions.INFILE_INDEX);
                    }
                } else {
                    ArrayList<String> indexFiles = new ArrayList<String>();
                    for(String sampleFile: props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP, 0)){
                        if(FileValidation.fileValid(sampleFile+".bai")){
                            indexFiles.add(sampleFile+".bai");
                        } else if(FileValidation.fileValid((sampleFile.replace(".bam", ".bai")))){
                            indexFiles.add(sampleFile.replace(".bam", ".bai"));
                        }else{
                            noBamIndex = true;
                            break;
                        }
                    }
                    if(noBamIndex){
                        logger.warn("Not all input files have a readable index file available. VAF determination will be disabled.");
                        props.remove(MobileDefinitions.INFILE_INDEX);
                    } else{
                        props.put(MobileDefinitions.INFILE_INDEX, String.join(",", indexFiles));
                    }
                }

                //When a resources directory has been provided, either in the properties file or on the command line,
                //it is checked whether it actually exists and an optional '/' is appended.
                if(props.containsKey(MobileDefinitions.RESOURCES_DIR)){
                    resourcesDir = props.getProperty(MobileDefinitions.RESOURCES_DIR);

                    //Check if exists
                    File dir = new File(resourcesDir);
                    if(!dir.isDirectory()){
                        printUsage();
                        logger.error("The provided resources folder does not exist.");
                        System.exit(-1);
                    }

                    //Add '/'
                    if(!resourcesDir.endsWith("/") && !resourcesDir.matches("(.|~)") && !resourcesDir.equals("")){
                        resourcesDir = resourcesDir + "/";
                        props.put(MobileDefinitions.RESOURCES_DIR, resourcesDir);
                    }

                } //Otherwise resourcesDir is set to an empty string
                else {
                    resourcesDir = "";
                    props.put(MobileDefinitions.RESOURCES_DIR, resourcesDir);
                }

                //Check if the reference genome is indeed available and in the correct format. If not, sequence determination will be disabled.
                if(props.containsKey(MobileDefinitions.REFERENCE_GENOME_FILE)){
                    try{
                        new ReferenceGenome(props.getProperty(MobileDefinitions.REFERENCE_GENOME_FILE));
                    } catch(IOException e){
                        logger.warn(e.getMessage() + ". Determination of the sequence for the target site duplication and insertion position will be disabled.");
                    }
                }

                //Check if predictions of a normal sample are provided. If not, filtering by normal predictions will be disabled.
                if(props.containsKey(MobileDefinitions.NORMAL_PREDICTIONS) && !FileValidation.fileValid(props.getProperty(MobileDefinitions.NORMAL_PREDICTIONS)))
                    logger.warn("The provided normal predictions file is not available or unreadable. Filtering by this file will be disabled.");
                else
                    logger.info("The provided normal predictions file will be used to filter the current somatic insertions.");

                //Check if 2 BAMs have been provided in case tumor insertion calling is set to true and if the name of the normal sample is present
                if(props.containsKey(MobileDefinitions.NORMAL_SAMPLE)){
                    String normalSampleName = props.getProperty(MobileDefinitions.NORMAL_SAMPLE);
                    if(sampleName.split(",").length != 2){
                        logger.warn("Less or more than 2 samples have been provided. Combined normal/tumor calling will be disabled.");
                        props.put(MobileDefinitions.NORMAL_TUMOR_CALLING, "false");
                    } else { props.put(MobileDefinitions.NORMAL_TUMOR_CALLING, "true"); }
                    if(!Arrays.asList(sampleName.split(",")).contains(normalSampleName)){
                        logger.warn("The provided name for the normal sample is not present among the other sample names. Combined normal/tumor calling will be disabled.");
                        props.put(MobileDefinitions.NORMAL_TUMOR_CALLING, "false");
                    } else { props.put(MobileDefinitions.NORMAL_TUMOR_CALLING, "true"); }

                    if("true".equals(props.getProperty(MobileDefinitions.NORMAL_TUMOR_CALLING)))
                        logger.info("Two samples and normal sample name detected. Combined normal/tumor calling is enabled.");
                }
            } else {
                printUsage();
                logger.error("No arguments provided. Please try again.");
                System.exit(-1);
            }
        } catch(ArrayIndexOutOfBoundsException e){
            printUsage();
            logger.error("Invalid arguments provided. Please try again. Make sure the correct arguments are followed by a value.");
            System.exit(-1);
        } catch (FileNotFoundException e) {
            logger.error(e.getMessage());
            System.exit(-1);
        } catch (IOException e) {
            logger.error(e.getMessage());
            System.exit(-1);
        }
    }

}