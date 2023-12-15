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

if [ $# -ne 1 ]; then echo -e "Script to extract coverage threshold.\nUsage: $0 <CSV>.\nVersion: 1.0"; exit 1; fi

export SUMMARY=$1

line=$(less ${SUMMARY} | head -n 1 | sed 's/,/\n/g' | awk '{print NR"\t"$1}' | grep mean_depth | cut -f1)
cov=$(less ${SUMMARY} | tail -n 1 | sed 's/,/\n/g' | head -n ${line} | tail -n1 | sed 's/\..*//g')
awk -v cov="$cov" 'BEGIN {print cov*12}'
