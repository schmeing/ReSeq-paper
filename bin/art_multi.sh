#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

cov=0
threads=48
output=""
params=""

while getopts "f:hj:o:p:" opt; do
  case "$opt" in
  f) cov=$OPTARG
    ;;
  h)
    echo "Usage: bash art_multi -j [threads] -f [total coverage] -o [output prefix] -p \"[art parameters]\""
    echo "Runs multiple intances of art and combines the results"
    exit 0
    ;;
  j) threads=$OPTARG
    ;;
  o) output=$OPTARG
    ;;
  p) params=$OPTARG
    ;;
  esac
done

shift $((OPTIND-1))
[ "$1" = "--" ] && shift

seq ${threads} | xargs -n1 -P${threads} -I@ bash -c "art_illumina ${params} -f $(echo "scale=6; $cov/$threads" | bc) -o ${output}_tmp@_ -d art@_ -na -rs \$RANDOM 1>${output}_tmp@.log 2>&1"
seq 2 | xargs -n1 -P2 -I@ bash -c "cat ${output}_tmp*_@.fq > ${output}@.fq"
rm -f ${output}_tmp*_?.fq 
