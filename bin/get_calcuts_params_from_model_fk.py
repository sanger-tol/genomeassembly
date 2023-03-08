#!/usr/bin/env python

import sys

def extract_cov(path):
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if line and line[0] == 'kmercov':
                return float(line[1])
    return -1
            

if __name__ == '__main__':
    path = sys.argv[1]
    cov = extract_cov(path)
    cov = int(cov + cov/2)
    max_val = cov*4
    calcuts_opts = '-m {:d} -u {:d}'.format(cov, max_val)
    print(calcuts_opts)
