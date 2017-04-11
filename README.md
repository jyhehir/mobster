+[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/mobster/README.html)
 +

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
java -Xmx8G -jar MobileInsertions-0.2.2.jar -properties Mobster.properties -out TestSample_

