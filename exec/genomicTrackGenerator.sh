#!/bin/bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# parser

nl=$'\n'

usage() { echo \
"$nl ##########  $(basename $0)  ########## $nl $nl usage: $nl $nl ex: genomicTrackGenerator.sh -r GATC -p 0 -n dpnII -g hg19.fa -o . $nl" && grep " .)\ #" $0;\
exit 0; }

[ $# -eq 0 ] && usage

while getopts ":hr:p:n:g:o:" arg; do
  case $arg in
    p) # restiction enzyme cut. For example pos 5 in HgaI GACGC indicates cleavage as follows: 5' GACGC^NNNNN
	pos=${OPTARG}
	;;
    r) # Restriction enzyme sequence
	re=${OPTARG}
	;;
    n) # Name restriction enzyme
	nameRe=${OPTARG}
	;;
    g) # Reference genome used for the generation of the genomic tracks
	refgen=${OPTARG}
	;;
    o) # Output path where to save the genomic tracks
	output=${OPTARG}
	;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

# set variables

nameRefgen=$(echo $(basename $refgen .fa))
genomicTracks=${output}/genomic_tracks_${nameRefgen}
splitChrom=${output}/splitChrom

# reorder pos array just in case 

# cut the restriction enzyme according to positions

cutRe=$re
cutRe=$(echo "${cutRe:0:${pos}}\t${cutRe:${pos}}")


# print variables 

echo Restriction enzyme  = "${re}"
echo Restriction enzyme name  = "${nameRe}"
echo Restiction enzyme cut  = "${pos}"
echo Genome will be cut by: $(echo $cutRe | sed 's/\\t/^/g')
echo Reference Genome = "${refgen}"
echo Output path = "${output}"

# resolve ambiguity on re

# standard abbreviations to represent ambiguity
# R = G or A
# Y = C or T
# M = A or C
# K = G or T
# S = G or C
# W = A or T
# B = not A (C or G or T)
# D = not C (A or G or T)
# H = not G (A or C or T)
# V = not T (A or C or G)
# N = A or C or G or T

cutRe=$(echo $cutRe | \
sed '{
	s/R/[GA]/g
	s/Y/[CT]/g
	s/M/[AC]/g
	s/K/[GT]/g
	s/S/[GC]/g
	s/W/[AT]/g
	s/B/[CGT]/g
	s/D/[AGT]/g
	s/H/[ACT]/g
	s/V/[ACG]/g
	s/N/[ACGT]/g
}')

cutReLow=$(echo $(echo $cutRe | tr '[:upper:]' '[:lower:]'))

# generate path

mkdir -p $genomicTracks
mkdir -p $splitChrom

# split genome by chrom

echo "`date` Generating genomic tracks of ${nameRefgen} using ${nameRe}..."

for i in $(grep ">" $refgen | sed "s/>//g");

do
	samtools faidx $refgen $i > $splitChrom/${i}.fa
done


# restriction enzyme to lowercase to find all restriction sites

reLow=$(echo "$re" | awk '{print tolower($0)}')


for file in $splitChrom/*.fa;
do
	chrL=$(cat $file | sed '1d' | tr -d "\n" | awk '{print length()}')
	chr=$(echo $(basename $file) | sed 's/.fa//')
	
	# define places where the r. site is present, correct the position and
	# sum the position where the enzyme cuts
	# generates the position for every chr separately
	
	echo "`date` Generating genomic tracks of ${chr} using ${nameRe}..."

	paste \
	<(cat $file | sed '1d' | tr -d "\n" | grep -io $reLow -aob | awk -v pos=$pos 'BEGIN {print "1"} {print $1+1+pos}') \
	<(cat $file | sed '1d' | tr -d "\n" | grep -io $reLow -aob | \
	awk -v pos=$pos -v chrL=$chrL '{print $1+pos};END{print chrL}') | \
	awk -v chr=$chr 'BEGIN {OFS="\t"}{print chr, $0}' >> \
	${genomicTracks}/${chr}_genomicTrackChr
	echo "`date` Finished generation of genomic tracks of ${chr} using ${nameRe}..."

done

# join and sort positions

cat $genomicTracks/*_genomicTrackChr | sort -k1,1 -k2,2n -k3,3n >> $genomicTracks/${nameRe}_genomicTrack

rm $genomicTracks/*_genomicTrackChr

echo "`date` Finished generation of genomic tracks of ${nameRefgen} using ${nameRe}, saved in ${genomicTracks}..."


