#!/bin/bash

FAI=$1

cut -f1-2 $FAI > $(basename $FAI).chrom.sizes
