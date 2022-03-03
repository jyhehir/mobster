package org.umcn.me.pairedend;

import java.io.*;
import java.util.Properties;
import java.util.Vector;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.umcn.me.output.Save;
import org.umcn.me.samexternal.SAMDefinitions;
import org.umcn.me.util.MobileDefinitions;

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

            //Filter and save the predictions
            Save.runFromProperties(props, predictions);

        } catch (IOException e) {
            logger.error(e.getMessage());
            System.exit(-1);
        }
    }

    // note that this is conditional, i.e. only if user requires this
    private static void configureAndRunPicard(final Properties props) throws IOException {
        //NOTE if multiple BAM Files are provided, only the insert size is investigated of the 1st BAM
        String picardCommand = "java -Xmx4g -jar " + resourcesDir + props.getProperty(MobileDefinitions.PICARD_COLLECT_INSERT_METRICS) +
                " VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_hist.pdf" +
                " INPUT=" + props.getProperty(MobileDefinitions.INFILE).split(MobileDefinitions.DEFAULT_SEP,0)[0] + " OUTPUT=" + props.getProperty(MobileDefinitions.OUTFILE).trim() + "_insertstats" +
                " STOP_AFTER=50000000";

        execUnixCommand(picardCommand);

        //Values are on line 8
        //[0 = MEDIAN] [4 = MEAN] [17=99percentile]
        BufferedReader br = new BufferedReader(new FileReader(props.getProperty(MobileDefinitions.OUTFILE) + "_insertstats"));
        String line;
        boolean read = false;
        double median;
        int mean = 0;
        double percentile99 = 0;
        int clustermax = 0;
        int sd = 0;

        while ((line = br.readLine()) != null){
            if(read){
                String[] split = line.split("\t", -1);
                median = Double.parseDouble(split[0]);
                mean = (int) Double.parseDouble(split[4]);
                percentile99 = Double.parseDouble(split[17]);
                clustermax = (int) (median + percentile99);
                sd = (int) Double.parseDouble(split[5]);

                break;
            }else if(line.startsWith("MEDIAN_INSERT")){
                read = true;
            }
        }

        br.close();

        if (mean == 0 || clustermax == 0 || sd == 0){
            logger.error("Could not parse the PICARD CollectInsertSizeMetrics successfully");
            System.exit(-1);
        }

        if (props.containsKey(MobileDefinitions.USE_READ_LENGTH) &&
                "true".equalsIgnoreCase(props.getProperty(MobileDefinitions.USE_READ_LENGTH).trim())){
            clustermax = clustermax - Integer.parseInt(props.getProperty(MobileDefinitions.READ_LENGTH).trim());
        }

        props.put(MobileDefinitions.LENGTH_99PROCENT_OF_FRAGMENTS, Integer.toString(clustermax));
        props.put(MobileDefinitions.MEAN_FRAGMENT_LENGTH, Integer.toString(mean));
        props.put(MobileDefinitions.SD_FRAGMENT_LENGTH, Integer.toString(sd));
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
        System.out.println("Only the properties file is required, but it is advisable to also provide the other single dash ('-[argument]') arguments.");
        System.out.println("These other arguments will override the default properties in the properties file.");
        System.out.println("In addition, any property as listed in the default properties file can also be overridden by providing it as a double dash argument ('--USE_PICARD [value]').");
        System.out.println("A new properties file containing all the final properties will also be created with the output prefix.");
        System.out.println();
        System.out.println("\t-properties [properties]\tThis value is required. Path to the properties file.");
        System.out.println(("\t-rdir [resources folder]\tThis value will override corresponding value in the properties file. Specify where the resources folder of mobster is installed."));
        System.out.println("\t-in [input .bam file]\t\tThis value will override corresponding value in the properties file. Multiple BAM files may be specified if separated by a comma");
        System.out.println("\t-out [output prefix]\t\tThis value will override corresponding value in the properties file.");
        System.out.println(("\t-sn [sample name]\t\t\tThis value will override corresponding value in the properties file. Multiple sample names may be specified if separated by a comma"));
        System.out.println(("\t-vcf\t\t\t\t\t\tIf this flag is given, the corresponding value in the properties file will be set to 'true'. Output the results of Mobster in VCF format in addition to the default format."));
        System.out.println(("\t--<PROPERTY> [value]\t\tThis value will override the corresponding value in the properties file."));
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
                        props.put(MobileDefinitions.VCF, true);
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
                    else if (argument.equalsIgnoreCase("-vcf")){
                        props.put(MobileDefinitions.VCF, value);
                    }
                    else if (argument.equalsIgnoreCase("-rdir")){
                        props.put(MobileDefinitions.RESOURCES_DIR, value);
                    }
                    //When a property present in the default properties is provided, it is overridden.
                    else if(argument.startsWith("--")
                            && props.containsKey(argUpperCase)){
                        props.put(argUpperCase, value);
                    } else if(!argument.equals("-properties")){
                        System.out.println(argument);

                        printUsage();
                        logger.error("Invalid arguments. Please try again.");
                        System.exit(-1);
                    }
                }

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

                //When a resources directory has been provided, either in the properties file or on the command line,
                //it is checked whether it actually exists and an optional '/' is appended.
                //In both cases, the final properties are written to a file
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

                    //Write the current properties to a file
                    String output = props.getProperty(MobileDefinitions.OUTFILE) + "_Mobster.properties";
                    props.store(new FileOutputStream(output), "");

                } //Otherwise resourcesDir is set to an empty string
                else {
                    //Write the current properties to a file
                    String output = props.getProperty(MobileDefinitions.OUTFILE) + "_Mobster.properties";
                    props.store(new FileOutputStream(output), "");

                    //Do not write the empty resourcesDir to a properties file
                    resourcesDir = "";
                    props.put(MobileDefinitions.RESOURCES_DIR, resourcesDir);
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