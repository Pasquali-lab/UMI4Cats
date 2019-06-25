#!/bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# parser

nl=$'\n'

usage() { echo \
"$nl ##########  $(basename $0)  ########## $nl $nl usage: $nl $nl ex: aligment.sh -w /imppc/labs/lplab/share/marc/umi4cBin/examples/prove/2019-06-18_umi4c_GTAC -i path/hg19.fa .b path/bowtie2  $nl" && grep " .)\ #" $0;\
exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hw:c:i:b:" arg; do
  case $arg in
    w) # Working directory for the UMI-4c analysis 
	wk_dir=${OPTARG}
	;;
    c) # Number of cores to use in the analysis
	cores=${OPTARG}
	;;
    i) # Reference genome in fasta format
	refgen=${OPTARG}
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

######### ALIGN ######### 

echo "`date` Starting aligment of files save into $wk_dir"

# set variables
out_alignment=${wk_dir}/alignment
index=$(echo ${refgen%.*})
mkdir -p $out_alignment

# alignment bowtie2
for file in ${wk_dir}/split/*_splited_clean.fastq;

do
	echo "`date` Aligning ${file}..."
	nonExt=$(basename $file _splited_clean.fastq)
	$bowtie2 -p $cores -x $index --reorder -U $file -S \
	${out_alignment}/${nonExt}.sam 2> ${out_alignment}/${nonExt}.log
	echo "`date` Finished alignment of ${file}..."	
done

echo "`date` Finished all aligments of files save into $wk_dir"








