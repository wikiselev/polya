#!/bin/bash
declare -a arr=(a-seq atlas)

rm out/*

# to check the jobs run: ps ax | grep rmd

for i in ${arr[@]}
do
	nohup /homes/kiselev/R-2.15.2/bin/Rscript rmd.R $i css-classification > out/$i-cc.out 2> out/$i-cc.err < /dev/null &
	nohup /homes/kiselev/R-2.15.2/bin/Rscript rmd.R $i protein-profiles > out/$i-pp.out 2> out/$i-pp.err < /dev/null &
done
