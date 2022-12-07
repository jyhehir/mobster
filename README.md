[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/mobster/README.html)

# mobster
NGS tool for detecting MEI and gene retrotransposition events in WGS and WES data, see Thung et al. Genome Biol. 2014 for more information.

## Install

### Dependencies

* Mobster depends on an aligner for aligning potential reads to mobiome reference. By default MOSAIK is the aligner of choice, which can be installed via bioconda or cloning the MOSAIK [git repository](https://github.com/wanpinglee/MOSAIK).

* Mobster is written in Java (tested working with Java 8), and built using maven, hence assumes that they are in your `PATH`.  

* `git lfs` is required. See [this](https://help.github.com/articles/installing-git-large-file-storage/) link for more detail.

### Simple installation

To build mobster simply run "install.sh". This will package the required classes into a fat executable jar in the "target" directory. Typically called MobileInsertions-\<version\>.jar. For example MobileInsertions-0.2.4.1.jar.

After installation, the aligner of choice (e.g. MOSAIK) is assumed to be in the `PATH`. If not, please don't forget to update Mobster.properties to include the location for the aligner.

## Testing the installation

Try running Mobster with the sample bam provided.

```bash
cd target
java -Xmx8G -jar MobileInsertions-0.2.5.0.jar \
    -properties Mobster.properties \
    -vcf \
    -out TestSample
```

Two relevant output files will be produced:
A table with predictions (.txt) and a VCF file.

## Single sample

Example on how to call the program for a single sample:

```bash
java -Xmx8G -jar MobileInsertions-0.2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in input.bam \
    -sn test_sample \
    -vcf \
    -out mobster_test
```

## Multiple samples

You can also run the program in multiple sample mode. For this you need to change `MULTIPLE_SAMPLE_CALLING` in the properties file to true. Then you can run Mobster like:

```bash
java -Xmx8G -jar MobileInsertions-2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in A1_child.bam,A1_father.bam,A1_mother.bam \
    -sn A1_child,A1_father,A1_mother \
    -vcf \
    -out A1_trio_mobster
```

Alternatively you can overrun any property in the property file by prepending any property name by '--' followed by the value:

```bash
java -Xmx8G -jar MobileInsertions-2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in A1_child.bam,A1_father.bam,A1_mother.bam \
    -sn A1_child,A1_father,A1_mother \
    -vcf \
    -out A1_trio_mobster \
    --MULTIPLE_SAMPLE_CALLING true
```

## Tumor/normal paired calling

Calling of somatic insertions in a paired tumor/normal sample is also possible through two different ways.
Either by providing both samples at the same time as when using multisample mode
(Additionally, the normal sample must be indicated with the '-ns' argument).

```bash
java -Xmx8G -jar MobileInsertions-2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in tumor.bam,normal.bam \
    -sn tumor,normal \
    -ns normal
    -vcf \
    -out tumor_normal
```

Another way is by first creating predictions on the normal sample

```bash
java -Xmx8G -jar MobileInsertions-2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in normal.bam
    -sn normal \
    -vcf \
    -out normal
```

Followed by providing those predictions with '-np' when running on the tumor sample:

```bash
java -Xmx8G -jar MobileInsertions-2.5.0.jar \
    -rdir ../resources
    -properties Mobster.properties \
    -in tumor.bam
    -sn tumor \
    -np normal_predictions.vcf
    -vcf \
    -out tumor_normal
```

Using the first method, the resulting insertion predictions will be labeled as either 'somatic', 'germline' or 'artifact' in the output txt file.
In the output VCF this is visible as either a PASS, germline or normal_artifact value respectively in the FILTER column.
The 'artifact' label is applied if the number of supporting reads in the normal sample is less than 3 as this could indicate possible noise or contamination of the tumor sample in the normal sample.
This threshold can be adjusted by changing the `MINIMUM_NORMAL_SUPPORTING_READS` property. 
However, as insertions in the normal sample are often missed in practice, any 'artifact' is often still an insertion in the normal sample. 
It is therefore often more practical to assume that only insertions labeled as 'somatic' are truly somatic.

Predictions made through the second method will only include predictions that are labeled as 'somatic'.


__Important__:

* Due to their size the test files and a number of resources for Mobster are stored with `git lfs`. These will not be retrieved with the initial git clone, but will be downloaded when you run "install.sh".
* If you want to turn on multiple sample calling, make sure you set the property `MULTIPLE_SAMPLE_CALLING` to "true" (the default is "false").
* You can set the `READ_LENGTH` property to the appropriate value for your BAM file in the properties file.
* In this version of Mobster, Mobster will accept BAM files produced by different mappers. In the property file `MAPPING_TOOL` is set to "unspecified" and `MINIMUM_MAPQ_ANCHOR` is set to 20. I would leave this as is.
* Most of the other default settings in the property file should be fine, to be a bit more sensitive you can lower the `MINIMUM_POLYA_LENGTH` to 7.