#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

threads=48
output=""
params=""
neat=""
samtools=""

while getopts "e:hj:o:p:s:" opt; do
  case "$opt" in
  e) neat=$OPTARG
    ;;
  h)
    echo "Usage: bash neat_multi.sh -j [threads] -c [total coverage] -e [neat folder] -o [output prefix] -p \"[neat parameters]\""
    echo "Runs multiple intances of neat and combines the results"
    exit 0
    ;;
  j) threads=$OPTARG
    ;;
  o) output=$OPTARG
    ;;
  p) params=$OPTARG
    ;;
  s) samtools=$OPTARG
    ;;
  esac
done

shift $((OPTIND-1))
[ "$1" = "--" ] && shift

rm -f ${output}_read{1,2}.fq ${output}_tmp_read{1,2}.fq.job*of*
seq ${threads} | xargs -n1 -P${threads} -I@ bash -c "python2 ${neat}/genReads.py ${params} -o ${output}_tmp --job @ ${threads} 1>${output}_tmp@.log 2>&1"
python2 ${neat}/mergeJobs.py -i ${output}_tmp -o ${output} -s ${samtools}

if test -f "${output}_read1.fq"; then
  rm -f ${output}_tmp_read{1,2}.fq.job*of*
else
  mv ${output}_tmp_read1.fq ${output}_read1.fq
  mv ${output}_tmp_read2.fq ${output}_read2.fq
fi
