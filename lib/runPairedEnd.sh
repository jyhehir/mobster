#!/bin/bash
#Paired End Detection of MEI Events
#NOTES:
# - network files of MOSAIK should be in same directory as shell script
# - jars should be in same directory as shell script
# - (shell script should be run from same directory?)
# - mosaik reference .dat and .jmp files should be in same directory as shell script
# - mosaik parameters can not (yet) be set through command line, set these by editing this script
#Author: Djie Thung

######Initialize variables##########
#When a numeric value has no default value it should be set to -1
fq1=""
fq2=""
in_bam=""
local_search=150
mean_fragment_length=-1
sd=-1
read_length=-1
hits_per_cluster=2
cluster_overlap=50
max_dist=700
search_area=-1
cluster_algorithm=""
out_prefix=""
samtools_dir=""
picard_dir=""
mosaik_dir=""
mobile_ref=""
human_ref=""
sep="/"
required_args=( "1" "2" "c" "o" "t" "g" ) #update when necessary
valid_algorithms=( "s" ) #update when necessary
########End of initializing variables############

function printUsage(){
	echo "############################################"
	echo "#Paired-End Mobile Element Insertion Finder#"
	echo "############################################"
	echo ""
	echo "Usage: pemeifinder [-12abcdfghiklmoprstz]"
	echo ""
	echo "OPTIONS:"
	echo "	-1 <fq file>		fq file for reads numbered 1 [REQUIRED]"
	echo "	-2 <fq file>		fq file for reads numbered 2 [REQUIRED]"
	echo "	-a <int>			search Area for simple clustering method [OPTIONAL, default = 2 * read length]"
	echo "	-b <BAM file>		bam file for calculating mean fragment length and sd [OPTIONAL]"
	echo "	-c <Cluster algorithm>		's' for simple clustering [REQUIRED]"
	echo "	-d <dir>		directory to samtools [OPTIONAL, default = '']"
	echo "	-f <int>		max distance between 5 and 3 end cluster for a double cluster (how Far) [OPTIONAL, default = $max_dist]"
	echo "	-g <dir>		path to human Genome reference archive and jump database (including file name stub) [REQUIRED]"
	echo "	-h <int>		min number of hits for a cluster [OPTIONAL, default = ${hits_per_cluster}]"
	echo "	-i <int>		max overlap (Intersection) for a 5 and 3 end cluster [OPTIONAL, default = $cluster_overlap]"
	echo "	-k <dir>		directory to mosaiK [OPTIONAL, default = '']"
	echo "	-l <int>		local search parameter. -l should be a number indicating how many times the standard deviation is searched. [OPTIONAL, default == $local_search]"
	echo "	-m <int>		mean fragment length (insert size) [OPTIONAL, but REQUIRED when -b not specified]"
	echo "	-o <prefix>		output prefix for all intermediate and final files. [REQUIRED]"
	echo "	-p <dir>		directory to Picard [OPTIONAL, default = '']"
	echo "	-r <int>		read length [OPTIONAL, default = $read_length]"
	echo "	-s <int>		standard deviation of fragment length distribution [OPTIONAL, but REQUIRED when -b not specified]"
	echo "	-t <dir>		path to transposable element reference archive and jump database (including file name stub) [REQUIRED]"
	echo "	-z <sep>		the Zeparator (seperator) used before pair numbers in args -1 and -2 [OPTIONAL, default = $sep]."
}

##Note: update this functions as arguments are added or removed (i know this is a lazy & dumb function right now!)
function printArguments(){
	echo ""
	echo "Arguments are set as follows: "
	echo "------------------------------"
	echo "fastq 1: ${fq1}"
	echo "fastq 2: ${fq2}"
	echo "original bam: ${in_bam}"
	echo "local search area: ${local_search}"
	echo "mean fragment length: ${mean_fragment_length}"
	echo "standard deviation of fragment length distribution: ${sd}"
	echo "read length: ${read_length}"
	echo "min hits per cluster: ${hits_per_cluster}"
	echo "max cluster overlap: ${cluster_overlap}"
	echo "max distance between 5 and 3 clusters: ${max_dist}"
	echo "max search area for simple clustering (may not be applicable): ${search_area}"
	echo "chosen cluster_algorithm: ${cluster_algorithm}"
	echo "output prefix: ${out_prefix}"
	echo "samtools dir (may be empty): ${samtools_dir}"
	echo "picard dir (may be empty): ${picard_dir}"
	echo "mosaik dir (may be empty): ${mosaik_dir}"
	echo "mobile ref file name: ${mobile_ref}"
	echo "human ref file name: ${human_ref}"
	echo "seperator for read groups: ${sep}"
}

#If param1 is not positive integer script will be exited
#param1: var to check whether it is a positive integer
#param2: name of var
function isInteger(){ 
	if ! [[ "$1" =~ ^[0-9]+$ ]]; then
		echo "error: value of -${2} should be positive integer" >&2
		exit 1
	fi
}

#If param1 is not positive integer or decimal script will be exited
#param1: var to check whether it is a positive integer or decimal
#param2: name of var
function isIntegerOrDecimal(){
	if ! [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
		echo "error: value of -${2} should be positive integer or decimal" >&2
		exit 1
	fi
}

if [ $# == 0 ]; then
	printUsage
	exit 0
fi

index=0
supplied_args={}

while getopts ":1:2:a:b:c:d:f:g:h:i:k:l:m:o:p:r:s:t:z:" opt; do
	supplied_args[index]=$opt
	let index=index+1
	case $opt in
		1)
			fq1=$OPTARG
			;;
		2)
			fq2=$OPTARG
			;;
		a)
			search_area=$OPTARG
			isInteger $search_area $opt
			;;
		b)
			in_bam=$OPTARG
			;;
		c)
			cluster_algorithm=$OPTARG
			if [ $(echo "${valid_algorithms[@]}" | grep -c $cluster_algorithm) == 0 ]; then
				echo "Invalid clustering algorithm specified." >&2
				echo "Valid arguments: ${valid_algorithms[@]}" >&2
				exit 1
			fi
			;;
		d)
			samtools_dir=$OPTARG
			;;
		f)
			max_dist=$OPTARG
			isInteger $max_dist $opt
			;;
		g)
			human_ref=$OPTARG
			;;
		h)
			hits_per_cluster=$OPTARG
			isInteger $hits_per_cluster $opt
			;;
		i)
			cluster_overlap=$OPTARG
			isInteger $cluster_overlap $opt
			;;
		k)
			mosaik_dir=$OPTARG
			;;
		l)
			#This should still be changed to local search = sd * -l
			local_search=$OPTARG
			isIntegerOrDecimal $local_search $opt	
			;;
		m)
			mean_fragment_length=$OPTARG
			isInteger $mean_fragment_length $opt
			;;
		o)
			out_prefix=$OPTARG
			;;
		p)
			picard_dir=$OPTARG
			;;
		r)
			read_length=$OPTARG
			isInteger $read_length $opt
			;;
		s)
			sd=$OPTARG
			isInteger $sd $opt
			;;
		t)
			mobile_ref=$OPTARG
			;;
		z)
			sep=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

###############################################
#Check whether all required arguments were set
###############################################
for i in "${required_args[@]}"; do
	if [ $(echo "${supplied_args[@]}" | grep -c $i) == 0 ]; then
		echo "-$i not supplied, required argument" >&2
		exit 1
	fi	
done

###################################################
#Check whether -m and -s are set when -b is not set
###################################################
if [ $(echo "${supplied_args[@]}" | grep -c b) == 0 ]; then
	if [ $(echo "${supplied_args[@]}" | grep -c m) == 0  ] ||
		[ $(echo "${supplied_args[@]}" | grep -c s) == 0 ]; then
		echo "-b not supplied, now -m and -s need to be supplied" >&2
		echo "set -m and -s " >&2
		exit 1
	fi
else
	java -Xmx4g -jar "${picard_dir}CollectInsertSizeMetrics.jar" VALIDATION_STRINGENCY=LENIENT TMP_DIR=./tmp HISTOGRAM_FILE="${out_prefix}_hist.pdf" INPUT="$in_bam" OUTPUT="${out_prefix}_insertstats"
	#Run picard here
	mean_fragment_length=$(sed -n '8p' "${out_prefix}_insertstats" | cut -f5)
	isIntegerOrDecimal $mean_fragment_length "m"
	mean_fragment_length=$(printf "%.0f" $(echo ${mean_fragment_length} | bc))
	sd=$(sed -n '8p' "${out_prefix}_insertstats" | cut -f6)
	isIntegerOrDecimal $sd "s"
	sd=$(printf "%.0f" $(echo ${sd} | bc))

fi	

if [ "$read_length" == -1 ]; then
	#!!!! IMPLEMENT CODE HERE TO ESTIMATE READ LENGTH
	echo "script for estimating read_length still needs to be developed" >&2
	echo "set -r for now, exiting..." >&2
	exit 1
fi

#If search_area is not set, set it to default of 2 * read_length
if [ "$search_area" == -1 ]; then
	let search_area=2*read_length
fi

printArguments

eval "${mosaik_dir}MosaikBuild -q ${fq1} -q2 ${fq2} -out ${out_prefix}.dat -st illumina -mfl ${mean_fragment_length} -quiet"
eval "${mosaik_dir}MosaikAligner -in ${out_prefix}.dat -out ${out_prefix} -ia ${human_ref}.dat -hs 15 -mm 10 -mhp 100 -act 35 -j ${human_ref} -p 4 -ls ${local_search} -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann -quiet"

java -Xmx20g -jar ./PotentialMEIFinder.jar -fq1 "$fq1" -fq2 "$fq2" -in "${out_prefix}.bam" -out "$out_prefix" -sep "$sep" 
cat "${out_prefix}_discordantunique.fq" "${out_prefix}_multiple.fq" "${out_prefix}_orphan.fq" > "${out_prefix}_potentialmobiles.fq"
eval "${mosaik_dir}MosaikBuild -q ${out_prefix}_potentialmobiles.fq -st illumina -out ${out_prefix}_potentialmobiles.dat -quiet"
eval "${mosaik_dir}MosaikAligner -in ${out_prefix}_potentialmobiles.dat -out ${out_prefix}_mappedpotentials -om -ia ${mobile_ref}.dat -hs 9 -mm 10 -act 20 -j ${mobile_ref}_hs9 -p 4 -annpe 2.1.26.pe.100.0065.ann -annse 2.1.26.se.100.005.ann -quiet"
eval "${samtools_dir}samtools sort -n ${out_prefix}_mappedpotentials.bam ${out_prefix}_mappedpotentials_nsorted"
eval "${samtools_dir}samtools sort -n ${out_prefix}_mappedpotentials.multiple.bam ${out_prefix}_mappedpotentials_multiple_nsorted"

java -Xmx20g -jar ./RefAndMEPairFinder.jar -multiple "${out_prefix}_mappedpotentials_multiple_nsorted.bam" -single "${out_prefix}_mappedpotentials_nsorted.bam" -sep "$sep" -potential "${out_prefix}_potential.bam" -out "$out_prefix"
eval "${samtools_dir}samtools sort ${out_prefix}_anchors.bam ${out_prefix}_anchors_csorted"
eval "${samtools_dir}samtools index ${out_prefix}_anchors_csorted.bam"
java -Xmx20g -jar ./AnchorClusterer.jar -cluster "$cluster_algorithm" -in "${out_prefix}_anchors_csorted.bam" -index "${out_prefix}_anchors_csorted.bam.bai" -maxdist "$max_dist" -mfl "$mean_fragment_length" -out "$out_prefix" -overlap "$cluster_overlap" -rl "$read_length" -rpc "$hits_per_cluster" -samdir "$samtools_dir" -sd "$sd" -search "$search_area"
