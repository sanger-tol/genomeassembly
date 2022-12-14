#!/bin/bash

args=$1
cpus=$2
buffer=$3
bed=$4
outfile=$5

sort $args --parallel=$cpus -S${buffer} $bed > $outfile
