#!/bin/bash

if [ $# -ne 1 ]; then echo -e "Script to extract coverage threshold.\nUsage: $0 <CSV>.\nVersion: 1.0"; exit 1; fi

export SUMMARY=$1

line=$(less ${SUMMARY} | head -n 1 | sed 's/,/\n/g' | awk '{print NR"\t"$1}' | grep mean_depth | cut -f1)
cov=$(less ${SUMMARY} | tail -n 1 | sed 's/,/\n/g' | head -n ${line} | tail -n1 | sed 's/\..*//g')
awk -v cov="$cov" 'BEGIN {print cov*12}'
