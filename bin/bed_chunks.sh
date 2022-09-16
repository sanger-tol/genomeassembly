#!/bin/bash

if [ $# -ne 2 ]; then echo -e "Script to split genome into chunks.\nUsage: $0 <FAI> <NUM_CHUNKS>.\nVersion: 1.0"; exit 1; fi

export FAI=$1
export CHUNKS=$2

for entry in $(awk '{print NR":"$1":"$2}' ${FAI}); do
    chunk_id=$(( $(echo $entry | cut -f 1 -d':') % ${CHUNKS}  ));
    contig=$(echo $entry | cut -f 2 -d':');
    end=$(echo $entry | cut -f 3 -d':');
    echo -e ${contig}"\t0\t"${end} >> ${chunk_id}.bed
done
