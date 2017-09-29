#!/bin/bash
#This shell script is a wrapper of IRA-SSMD.R for pairwise score calculation
#Due to the row number variation among distinct control-treatment pairs, one table could not process all the data at one time 
if [ $# -eq 0 ]; then echo "usage: $0 RScriptname[mandatory] table[mandatory] output_prefix[optional,default output] threshold[optional,default 10] indicator[to use sum of ctrl+sample(specify 1) for filtration, default 0(use ctrl only)]"; exit 1 ;
elif [ $# -eq 1 ]; then echo "R script is defined, plese specify the input"; exit 1; 
elif [ $# -eq 2 ]; then rscript=${1} bigTable=${2};outTag="output";threshold=10; sumfun=0;
elif [ $# -eq 3 ]; then rscript=${1} bigTable=${2};outTag=${3};threshold=10; sumfun=0; 
elif [ $# -eq 4 ]; then rscript=${1} bigTable=${2};outTag=${3};threshold=${4}; sumfun=0; 
else 
rscript=${1}; bigTable=${2};outTag=${3};threshold=${4};sumfun=${5};
fi
nameString=$(head -n 1 $bigTable| tr -d '\r')
IFS=', ' read -r -a namearray <<< $nameString
echo $nameString;
cols=$(awk 'BEGIN{FS=","}; END{print NF}' $bigTable)
treats=$(expr $cols - 3)
echo $treats treatments in the table...
echo "output prefix is $outTag"
echo "Read counts lower than $threshold will be removed before SSMD score calculation"
echo "get pairwise tables and calculate the SSMD score"
for i in $(seq 4 $cols); do cut -d , -f1,2,3,$i $bigTable > ${namearray[2]}_${namearray[$i-1]}.csv; done;
for i in $(seq 4 $cols); do Rscript --vanilla ${rscript} ${namearray[2]}_${namearray[$i-1]}.csv ${outTag}_${namearray[2]}_${namearray[$i-1]} $threshold ${sumfun}; done; 
