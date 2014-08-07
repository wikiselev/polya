#!/bin/bash

# declare -a arr=(element1 element2 element3)
declare -a arr=(a-seq atlas)

for i in ${arr[@]}
do
	bsub /homes/kiselev/R-2.15.2/bin/Rscript protein-profiles.R $i
done
