#!/bin/bash

export FAI=$1
export CHUNKS=$2

for entry in $(cat ${FAI} | awk '{print NR":"$1":"$2}' ); do
    chunk_id=$(echo $entry | cut -f 1 -d':');
    chunk_id=$(( $chunk_id % ${CHUNKS} ));
    echo ${GROUPS} $chunk_id
    contig=$(echo $entry | cut -f 2 -d':');
    end=$(echo $entry | cut -f 3 -d':');
    echo -e ${contig}"\t0\t"${end} >> ${chunk_id}.bed
done
