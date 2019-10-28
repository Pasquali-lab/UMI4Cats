#!/bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# parser

nl=$'\n'

usage() { echo \
"$nl ##########  $(basename $0)  ########## $nl $nl usage: $nl $nl ex: split2.sh -w /imppc/labs/lplab/share/marc/umi4cBin/examples/prove/2019-06-18_umi4c_GTAC -c 3 -r GTAC -t path/trimmomatic  $nl" && grep " .)\ #" $0;\
exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hw:c:r:t:" arg; do
  case $arg in
    w) # Working directory for the UMI-4c analysis 
	wk_dir=${OPTARG}
	;;
    c) # Number of cores to use in the analysis
	cores=${OPTARG}
	;;
    r) # Restriction enzyme sequence
	re=${OPTARG}
	;;
    t) # Path trimmomatic
	trimmomatic=${OPTARG}
	;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

######### SPLIT ######### 

echo "`date` Starting splitting of files save into $wk_dir"

# set variables

threads=$(echo "$(($cores * 2))")
mkdir -p ${wk_dir}/split

# split reads into fragments according to the restriction enzyme
# then qualitiy check

echo "`date` Spliting reads according to restriction enzyme ${re} and applying quality check..."

for file in ${wk_dir}/prep/*_filtered_R1.fastq;
do
	file_R1=$file
	file_R2=$(echo $file_R1 | sed 's/_R1./_R2./')
	nonExt=$(basename ${file%_*_*})
	outSplit=${wk_dir}/split/${nonExt}_splited.fastq
	outTrimo=${wk_dir}/split/${nonExt}_splited_clean.fastq

	echo "`date` Starting ${nonExt}..."

	# concat files, split every read by its restriction enzyme and generated a marged 
	# fastq with every new fragment
	# split qc line according to the length of the read
	# in the header is added: SEG (number of frag):starting pos of frag: end pos of frag, we count the lenght of re

	# generate the first column with the length of every read, starting with 1000 bp in R2
	# attach sequence and quality to the header
	# print as a fastq 

	paste \
	<(cat $file_R1 $file_R2 | paste - - - - | awk -v re=$re 'BEGIN {FS="\t"; OFS="\t"}\
	{gsub(re,"\t", $2); {split ($2, a, "\t")};\
	for (i in a) print $1":SEG"NR":"length(a[i-1])":"length(a[i-1])+length(a[i]),\
	a[i], $3, substr($4, 1, length(a[i]))}' |\
	awk 'BEGIN {FS="\t"; OFS="\t"} {print $1}'|\
	awk 'BEGIN {FS=":"; OFS=":"} {if ($9 == "R2") {$12 = 1000 - $12; $11 = 1000 - $11} print $0}') \
	<(cat $file_R1 $file_R2 | paste - - - - |  awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $3, $4}')|\
	sed 's/\t/\n/g' > $outSplit

	# clips the read once the average quality within the windowfalls below 
	# 20 and cut it, keep the once with minimum lenght of 30
	
	java -jar $trimmomatic SE -phred33 -threads $threads $outSplit $outTrimo SLIDINGWINDOW:4:20 MINLEN:20
	echo "`date` Finished ${nonExt}..."
done

