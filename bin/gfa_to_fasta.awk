## Convert GFA to FASTA format
## Date: 2025-04-04

BEGIN { OFS = "\t" }
/^S/ {
    print ">" $2
    print $3
}
