# mobster
NGS tool for detecting MEI and gene retrotransposition events in WGS and WES data


## Publication:

Mobster: accurate detection of mobile element insertions in next generation sequencing data.

[Thung et al. Genome Biol. 2014](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-014-0488-x.pdf)


## Install

### Dependencies

* Mobster depends on an aligner for aligning potential reads to mobiome reference. By default MOSAIK is the aligner of choice, which can be installed via bioconda or cloning the MOSAIK [git repository](https://github.com/wanpinglee/MOSAIK).

* Mobster is written in Java (tested working with Java 8), and built using maven, hence assumes that they are in your `PATH`.  

* `git lfs` is required. See [this](https://help.github.com/articles/installing-git-large-file-storage/) link for more detail.

### Simple installation

The sources can be cloned to any directory:

```
git clone git@github.com:jyhehir/mobster.git
```

Then, to build mobster simply run "install.sh". 

```
cd mobster
./install.sh
```

This will package the required classes into a fat executable jar in the "target" directory. Typically called MobileInsertions-\<version\>.jar. For example MobileInsertions-0.2.4.1.jar.

After installation, the aligner of choice (e.g. MOSAIK) is assumed to be in the `PATH`. If not, please don't forget to update Mobster.properties to include the location for the aligner.

For GRCh38 use, the user should unpack the compressed repmask resources:
```
gunzip alu_l1_herv_sva_other_grch38_accession_ucsc.rpmsk.gz
```


## Testing the installation

Try running Mobster with the sample bam provided

```bash
cd target
java -Xmx8G -jar MobileInsertions-0.2.4.1.jar \
    -properties Mobster.properties \
    -out TestSample
```

Example on how to call the program for a single sample:

```bash
java -Xmx8G -jar MobileInsertions-1.0-SNAPSHOT.jar \
    -properties Mobster_latest.properties \
    -in input.bam \
    -sn test_sample \
    -out mobster_test
```

You can also run the program in multiple sample mode. For this you need to change `MULTIPLE_SAMPLE_CALLING` in the properties file to true. Then you can run Mobster like:

```bash
java -Xmx8G -jar MobileInsertions-1.0-SNAPSHOT.jar \
    -properties Mobster_latest.properties \
    -in A1_child.bam,A1_father.bam,A1_mother.bam \
    -sn A1_child,A1_father,A1_mother \
    -out A1_trio_mobster
```

__Important__:

* Due to their size the test files and a number of resources for Mobster are stored with `git lfs`. These will not be retrieved with the initial git clone, but will be downloaded when you run "install.sh".
* If you want to turn on multiple sample calling, make sure you set the property `MULTIPLE_SAMPLE_CALLING` to "true" (the default is "false").
* You can set the `READ_LENGTH` property to the appropriate value for your BAM file in the properties file.
* In this version of Mobster, Mobster will accept BAM files produced by different mappers. In the property file `MAPPING_TOOL` is set to "unspecified" and `MINIMUM_MAPQ_ANCHOR` is set to 20. I would leave this as is.
* Most of the other default settings in the property file should be fine, to be a bit more sensitive you can lower the `MINIMUM_POLYA_LENGTH` to 7.

## RepeatMasker file

Mobster needs a GRCh37/GRCh38 repmask library file.

### GRCh37/GRCh38
The location of this file is included in the `Mobster.properties` file:
- GRCh37 (default):
`REPEATMASK_FILE=../resources/repmask/hg19_alul1svaerv.rpmsk`
- The user can change this default value to GRCh38:
`REPEATMASK_FILE=../resources/repmask/alu_l1_herv_sva_other_grch38_ucsc.rpmsk`

WARNING: User should unpack the compressed repmask resources (alu_l1_herv_sva_other_grch38_accession_ucsc.rpmsk.gz) during the install.

### Updating the RepeatMasker file (if needed)
You can freely download an updated repmask file from the "http://genome.ucsc.edu/cgi-bin/hgTables". 
There are many output options, here are the changes that you'll need to make:
- “GRCh37” or “GRCh38” assembly
- "Repeats" group
- "Repeatmasker" track
- “genome” region. 
- Select output format as "selected fields from primary and related tables"
- Choose the following output filename: `alu_l1_herv_sva_other_grch37_ucsc.rpmsk.tmp` or `alu_l1_herv_sva_other_grch38_ucsc.rpmsk.tmp`. 
- Cick the "get output" button
- Select all fields (click the "check all" button)
- Cick the "get output" button

Then, you have to filter this file and keep only the lines with:
- repFamily = 'Alu' OR 'L1' OR 'SVA' 
- or repName like 'HERV%' 
