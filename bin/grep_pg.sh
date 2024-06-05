#!/bin/bash

# grep_pg.sh
# -------------------
# A shell script to exclude pg lines and label read 1 and read 2 from cram containers
#
# -------------------
# Author = yy5

grep -v "^@PG" | awk '{if($1 ~ /^@/) {print($0)} else {if(and($2,64)>0) {print("1"$0)} else {print("2"$0)}}}'
