#!/usr/bin/env python
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

import sys


def extract_cov(path):
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if line and line[0] == "kmercov":
                return float(line[1])
    return -1


if __name__ == "__main__":
    path = sys.argv[1]
    cov = extract_cov(path)
    cov = int(cov + cov / 2)
    max_val = cov * 4
    calcuts_opts = "-m {:d} -u {:d}".format(cov, max_val)
    print(calcuts_opts)
