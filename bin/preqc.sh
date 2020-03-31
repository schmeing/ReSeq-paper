#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

output="preqc/preqc"
logfile="preqc/preqc.log"
threads=64

while getopts "h?j:l:o:" opt; do
  case "$opt" in
  h|\?)
    echo "Usage: bash preqc.sh [OPTIONS] File1 File2 ..."
    echo "Creates a preqc report for paired-end data in File1 and File2"
    echo "  -h -?                    display this help"
    echo "  -j [INT]                 defines the number of threads run"
    echo "  -l [STRING]              defines log file [preqc/preqc.log]"
    echo "  -o [STRING]              defines output prefix for report files [preqc/preqc]"
    exit 0
    ;;
  j) threads=$OPTARG
    ;;
  l) logfile=$OPTARG
    ;;
  o) output=$OPTARG
    ;;
  esac
done

shift $((OPTIND-1))
[ "$1" = "--" ] && shift

mkdir -p $(dirname "${output}")

echo "preprocess:" > ${logfile} &&\
sga preprocess --pe-mode 1 $@ > ${output}.fastq 2>>${logfile} &&\
\
echo "index:" >> ${logfile} &&\
sga index -a ropebwt --no-reverse -t ${threads} -p ${output} ${output}.fastq >>${logfile} 2>&1 &&\
\
echo "preqc:" >> ${logfile} &&\
sga preqc -t ${threads} ${output}.fastq > ${output}.preqc 2>>${logfile}

rm -f ${output}.fastq ${output}.bwt ${output}.sai
