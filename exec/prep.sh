#!/bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# parser

nl=$'\n'

usage() { echo \
"$nl ##########  $(basename $0)  ########## $nl $nl usage: $nl $nl ex: prep.sh -i path/raw -c 4 -s GTTGTCCTTGGGTTTAGCTGC -p ACCTC -r GTAC -g path/hg19.fa -w . -f fastq-multx $nl" && grep " .)\ #" $0;\
exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hi:c:s:p:r:g:w:f:" arg; do
  case $arg in
    i) # path with the fastq to analyse
	wk_path=${OPTARG}
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
    g) # Reference genome in fasta format
	refgen=${OPTARG}
	;;
    w) # Working directory for the UMI-4c analysis 
	wk_dir=${OPTARG}
	;;
    f) # Fastq-multx path
	fastqmultx=${OPTARG}
	;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

######### PREPARE ######### 

echo "`date` Starting preparation of files save into ${wk_path}"

# set variables

mkdir -p ${wk_dir}/prep

# Check for gz files and extract them

gz_files=$(find $wk_path -type f -name "*.gz")


if [ ! -z "$gz_files" ]
then
	echo "`date` Decompressing $gz_files"
	for file in $gz_files; do echo "gunzip -c $file > ${wk_dir}/prep/$(basename $file .gz)"; done | parallel

else 
	echo "`date` Non files to decompress"
fi

# generate umi identifier

for file in $(find $wk_path -maxdepth 1 -name "*_R2.fastq" -or -name "*_R2.fq");
do
	# R2 with umi in the header
	
	nonExt=$(basename ${file%_*})
	R2_umi=${wk_dir}/prep/${nonExt}_umi_R2.fastq
	cat $file | paste - - - - | \
	awk 'BEGIN {FS="\t"; OFS="\t"} {split($1, a, " ")}; {print a[1]":"substr($2, 1, 10), $2, $3, $4}' | \
	sed 's/\t/\n/g' > $R2_umi

	#R1 with umi in the header
	
	R1_umi=${wk_dir}/prep/${nonExt}_umi_R1.fastq
	paste \
	<(cat $R2_umi| paste - - - - | awk 'BEGIN {FS="\t"} {print $1}') \
	<(cat ${wk_path}/${nonExt}_R1.fastq | paste - - - - | awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $3, $4}') | \
	sed 's/\t/\n/g'	> $R1_umi

done

# Parse only the baits' reads and filter out non-specific reads

# create file with barcode using bait_seq, bait_pad and re 
# reads that present the bait and re will be keeped in the *_filtered* file

printf "prefiltered\t${bait_seq}${bait_pad}${re}" > ${wk_dir}/prep/barcodes.txt

echo "`date` Filtering reads..."


# loop taking on account both formats, just in case

for file in ${wk_dir}/prep/*_umi_R1.fastq;
do 
	file_R1=$file
	file_R2=$(echo $file_R1 | sed 's/_R1./_R2./')
	nonExt=$(basename ${file%_*_*})

	echo "`date` Starting ${nonExt}..."

	${fastqmultx} \
		-x \
		-m 0 \
		-b ${wk_dir}/prep/barcodes.txt $file_R1 $file_R2 \
		-o $wk_dir/prep/${nonExt}_%_R1.fastq \
		-o $wk_dir/prep/${nonExt}_%_R2.fastq
	
	outR1=$(echo $wk_dir/prep/${nonExt}_prefiltered_R1.fastq)
	outR2=$(echo $wk_dir/prep/${nonExt}_prefiltered_R2.fastq)
 	out2R1=$(echo $wk_dir/prep/${nonExt}_filtered_R1.fastq)
	out2R2=$(echo $wk_dir/prep/${nonExt}_filtered_R2.fastq)

	# remove annoying header modification done by fastq-multx
	# r1 and r2 is added for each file respectively
	cat $outR1 | paste - - - - | \
	awk 'BEGIN {FS="\t"; OFS="\t"}; {split($0, a, " ")}; {print a[1]":R1", $2, $3, $4}' | \
	sed 's/\t/\n/g' > $out2R1

	cat $outR2 | paste - - - - | \
	awk 'BEGIN {FS="\t"; OFS="\t"}; {split($0, a, " ")}; {print a[1]":R2", $2, $3, $4}' |
	sed 's/\t/\n/g' > $out2R2

	echo "`date` Finished ${nonExt}..."
done

echo "`date` Finished preparation of files save into ${wk_path}"


