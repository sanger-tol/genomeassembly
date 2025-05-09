## Script to split genome into chunks
##
## Author: Jim Downie

BEGIN { OFS = "\t" }
{
    chunk_id = NR % chunk_size
    print $1, 0, $2 > prefix "." chunk_id ".bed"
}
