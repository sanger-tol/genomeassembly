#!/bin/bash

#
# Based on https://github.com/sanger-tol/treeval/blob/80554a803903183613d49690d5770eeadb3c42c9/bin/generate_cram_csv.sh
# from Sanger TOL treeval pipeline
#

chunkn=0
for cram in "$@"; do

    rgline=$(samtools view -H $cram|grep "RG"|sed 's/\t/\\t/g'|sed "s/'//g")

    crampath=$(readlink -f ${cram})

    ncontainers=$(zcat ${crampath}.crai|wc -l)
    base=$(basename $cram .cram)

    from=0
    to=10000


    while [ $to -lt $ncontainers ]
    do
        echo $crampath,${crampath}.crai,${from},${to},${base},${chunkn},${rgline}
        from=$((to+1))
        ((to+=10000))
        ((chunkn++))
    done

    if [ $from -le $ncontainers ]
    then
        echo $crampath,${crampath}.crai,${from},${ncontainers},${base},${chunkn},${rgline}
        ((chunkn++))
    fi
done

exit 0
