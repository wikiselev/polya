#!/bin/bash

# declare -a arr=(element1 element2 element3)
declare -a arr=(a-seq atlas)

for i in ${arr[@]}
do
	bsub -J $i /homes/kiselev/R-2.15.2/bin/Rscript $i.R
	# remove old mappings
	rm -r ../../data/map/$i
	mkdir ../../data/map/$i
	# count number of processed protein files
	targetDir=../../data/prot/processed/
	n=`find ${targetDir} -type f | wc -l`
	# run a array job for each protein file
	bsub -J "$i-map[1-$n]" -w "ended("$i")" /homes/kiselev/R-2.15.2/bin/Rscript mapping.R $i
	# merge all mappings in one file
	bsub -w "ended("$i-map")" /homes/kiselev/R-2.15.2/bin/Rscript merge_mapping.R $i
done
