#!/usr/bin/env bash
echo "Downloading additional resources"
git lfs pull
gunzip ./resources/repmask/alu_l1_herv_sva_other_grch38_ucsc_repeatmasker.txt.gz

echo -e "${GREEN} downloading human mobiome files${NC}"

mkdir mobiome
wget -P mobiome https://github.com/jyhehir/mobster/raw/master/resources/mobiome/54_mobiles_inclHERVK.dat
wget -P mobiome https://github.com/jyhehir/mobster/raw/master/resources/mobiome/54_mobiles_inclHERVK_hs9_keys.jmp
wget -P mobiome https://github.com/jyhehir/mobster/raw/master/resources/mobiome/54_mobiles_inclHERVK_hs9_meta.jmp
wget -P mobiome https://github.com/jyhehir/mobster/raw/master/resources/mobiome/54_mobiles_inclHERVK_hs9_positions.jmp

echo -e "${GREEN}downloading MOSAIK network files${NC}"
mkdir MOSAIK
wget -P MOSAIK https://github.com/jyhehir/mobster/raw/master/resources/MOSAIK/2.1.26.pe.100.0065.ann
wget -P MOSAIK https://github.com/jyhehir/mobster/raw/master/resources/MOSAIK/2.1.26.se.100.005.ann

echo -e "${GREEN}downloading picard-1.73 dependency${NC}"
mkdir picard-1.73
wget -P picard-1.73 https://github.com/jyhehir/mobster/raw/master/resources/picard-1.73/CollectInsertSizeMetrics.jar

echo -e "${GREEN}downloading repeatmasker files${NC}"
mkdir repmask
wget -P repmask https://github.com/jyhehir/mobster/raw/master/resources/repmask/hg19_alul1svaerv.rpmsk
mv repmask/hg19_alul1svaerv.rpmsk ${script_dir}/repmask/hg19_alul1svaerv.txt
wget -P repmask https://github.com/jyhehir/mobster/raw/develop/resources/repmask/alu_l1_herv_sva_other_grch38_ucsc_repeatmasker.txt.gz

gunzip repmask/alu_l1_herv_sva_other_grch38_ucsc_repeatmasker.txt.gz

echo -e "${GREEN}downloading test input files${NC}"
mkdir test_files

wget -P test_files https://github.com/jyhehir/mobster/raw/develop/test_files/chr12_wgs_3000ins_100excl_10x_bwadef.bam
mkdir test_output

#MOSAIK INSTALLATION
echo -e "${GREEN}Downloading of Mobster dependencies and resources has been completed${NC}\n\n"

echo -e "${RED}NOTE${NC}: You will still have to manually install MOSAIK version ${RED}2.1.33${NC} and add the binaries to your ${RED}PATH${NC}\n"
echo -e "You can download the x64 binaries or the source from: https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/mosaik-aligner/MOSAIK-2.1.33-Linux-x64.tar or https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/mosaik-aligner/MOSAIK-2.1.33-source.tar\n"

echo "Starting mobster build..."

mvn clean

cd ./lib/jar

mvn install:install-file -Dfile=java-yield-1.0-SLIMJAR.jar -DgroupId=java-yield -DartifactId=java-yield -Dversion=1.0-SLIMJAR -Dpackaging=jar
mvn install:install-file -Dfile=picard-1.113.jar -DgroupId=net.sf.picard -DartifactId=picard -Dversion=1.113 -Dpackaging=jar
mvn install:install-file -Dfile=samtools-1.113.jar -DgroupId=net.sf.samtools -DartifactId=samtools -Dversion=1.113 -Dpackaging=jar

cd ../..

mvn compile -U
mvn package

cp ./lib/Mobster.properties ./target
