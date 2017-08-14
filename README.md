[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/mobster/README.html)

# mobster
NGS tool for detecting MEI and gene retrotransposition events in WGS and WES data, see Thung et al. Genome Biol. 2014 for more information.

Mobster is built using maven, and assumes that this is in your path. To build mobster run "install.sh" this will package the required classes into a fat executable jar in the "target" directory. Typically called MobileInsertions-.jar. For example MobileInsertions-0.2.2.jar.

#Requirements
Mobster also requires MOSAIK to be installed,
this can be done via
bioconda or git clone https://github.com/wanpinglee/MOSAIK.git

Don't forget to update Mobster.properties to include the location for MOSAIK

#Testing the install

Try running Mobster with the sample bam provided

java -Xmx8G -jar MobileInsertions-0.2.2.jar -properties Mobster.properties -out TestSample

Important:
* Its best to use the properties file Mobster_latest.properties contained within the zip of the modified version of Mobster.
* if you want to turn on multiple sample calling, make sure you set the property MULTIPLE_SAMPLE_CALLING=false to true. (the default is false)
* You can set the READ_LENGTH property to whatever value is appropriate for your BAM file in the properties file
* In this version of Mobster, Mobster will accept all kinds of BAM files from different mappers. In the property file MAPPING_TOOL is set to unspecified and MINIMUM_MAPQ_ANCHOR is set to 20. I would leave this as is.
* Most of the other default settings in the property file should be fine, to be a bit more sensitive you can lower the MINIMUM_POLYA_LENGTH to 7

Example on how to call the program for a single sample:

java -Xmx8G -jar MobileInsertions-1.0-SNAPSHOT.jar -properties Mobster_latest.properties -in input.bam -sn test_sample -out mobster_test

You can also run the program in multiple sample mode. For this you need to change MULTIPLE_SAMPLE_CALLING in the properties file to true. Then you can run Mobster like:

java -Xmx8G -jar MobileInsertions-1.0-SNAPSHOT.jar -properties Mobster_latest.properties -in A1_child.bam,A1_father.bam,A1_mother.bam -sn A1_child,A1_father,A1_mother -out A1_trio_mobster
