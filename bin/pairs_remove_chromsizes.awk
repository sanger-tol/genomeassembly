## Remove the #chromsize: lines from
## a pairs file, to make it compatible with
## JUICERTOOLS_PRE.
##
## Author: Jim Downie

BEGIN { FS = OFS = "\t" }
!/^#chromsize/ {
    print $0
}
