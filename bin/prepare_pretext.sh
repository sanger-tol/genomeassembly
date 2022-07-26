#!/bin/bash

FAI=$1
PAIRS=$2
awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END{print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' $FAI ; 
awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' $PAIRS
