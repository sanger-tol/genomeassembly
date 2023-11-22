#!/bin/bash
#
#    Copyright (C) 2022-2023 Genome Research Ltd.
#
#    Author: Ksenia Krasheninnikova <kk16@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

if [ $# -ne 2 ]; then echo -e "Script to split genome into chunks.\nUsage: $0 <FAI> <NUM_CHUNKS>.\nVersion: 1.0"; exit 1; fi

export FAI=$1
export CHUNKS=$2

for entry in $(awk '{print NR":"$1":"$2}' ${FAI}); do
    chunk_id=$(( $(echo $entry | cut -f 1 -d':') % ${CHUNKS}  ));
    contig=$(echo $entry | cut -f 2 -d':');
    end=$(echo $entry | cut -f 3 -d':');
    echo -e ${contig}"\t0\t"${end} >> ${chunk_id}.bed
done
