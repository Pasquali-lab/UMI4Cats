#!/bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# parser

nl=$'\n'

usage() { echo \
"$nl ##########  $(basename $0)  ########## $nl $nl usage: $nl $nl ex: umiCounter.sh -w 2019-06-19_umi4c_GTAC -c 3 -s GTTGTCCTTGGGTTTAGCTGC -p ACCTCT -r GTAC -i /hg19/hg19.fa -f genomic_tracks_hg19/csp6i_genomicTrack  -b path/bowtie2 $nl" && grep " .)\ #" $0;\
exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hw:c:s:p:r:i:f:b:" arg; do
  case $arg in
    w) # Working directory for the UMI-4c analysis 
	wk_dir=${OPTARG}
	;;
    c) # Number of cores to use in the analysis
	cores=${OPTARG}
	;;
    s) # Bait sequence
	bait_seq=${OPTARG}
	;;
    p) # Bait pad
	bait_pad=${OPTARG}
	;;
    r) # Restriction enzyme sequence
	re=${OPTARG}
	;;
    i) # Reference genome in fasta format
	refgen=${OPTARG}
	;;
    f) # Genomic track path
	genomicTrack=${OPTARG}
	;;
    b) # Path bowtie
	bowtie2=${OPTARG}
	;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

######### UMI COUNTER ######### 

echo "`date` Starting umi counter of files save into ${wk_dir}/alignment"

# set variables

index=$(echo ${refgen%.*})
mkdir -p $wk_dir/umiCounter
bowtie2=bowtie2

# set coordinates of the bait taking only perfect matches

coor=($($bowtie2 --quiet -x $index -c $(echo ${bait_seq}${bait_pad}${re}) -N 0 | \
samtools view | \
awk '{print $3,$4}'))

chr=$(echo ${coor[0]})
start=$(echo ${coor[1]})

for sam in ${wk_dir}/alignment/*.sam;
do
	fileNonExt=$(basename ${sam%.*})
	echo "`date` Starting with ${sam}..."
	# create a file that relay segments and fragments

	# select only segments aligned with more than 42 mapq, in the chr of bait and well aligned

	# intersec segments and frags counting the number of matched nucleotides
	bedtools intersect \
	-wo -a <(cat $genomicTrack | awk -v chr=$chr 'BEGIN {OFS="\t"} $1 == chr {print $1,$2,$3}') \
	-b <(samtools view $sam | awk -v chr=$chr 'BEGIN {OFS="\t"} (($2 == 0 || $2 == 16) && ($3 == chr) && ($5>=42)) {print $3,$4,$4+length($10)}') | uniq | \
	# sort and keep the frag-seg with more matched nucleotides
	sort -nk5.5 -nk6.6 -nrk7.7 | \
	awk '!_[$4,$5,$6]++ {print $4,$5,$6,$1,$2,$3}' > ${wk_dir}/umiCounter/${fileNonExt}SegFrag.txt


	# algorithm umi filtering???

	# count umis
	# frag-seg 
	# file with alignment and umi information and segment coor
	# collapse umis
	# recount umis

	join \
	<(awk '{print $1"-"$2"-"$3,$4,$5,$6}' ${wk_dir}/umiCounter/${fileNonExt}SegFrag.txt | sort) \
	<(samtools view $sam | awk -v chr=$chr 'BEGIN {OFS="\t"} (($2 == 0 || $2 == 16) && ($3 == chr) && ($5>=42)) {print $3"-"$4"-"$4+length($10),$1,$5,$10}' | sort) | \
	awk -v start=$start -v chr=$chr 'BEGIN {OFS="\t"} {split ($5,a,":")} {print chr, start, $2, $3, a[8]}' | uniq | cut -f1,2,3,4 | \
	uniq -c | sort -rnk1,1 | awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5,$1}' > ${wk_dir}/umiCounter/${fileNonExt}UmiRecount.txt

	echo "`date` Finished with ${sam}..."
done

echo "`date` Finished umi counter of files save into ${wk_dir}"
