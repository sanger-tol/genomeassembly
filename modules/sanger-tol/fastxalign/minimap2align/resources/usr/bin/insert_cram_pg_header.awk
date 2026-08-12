#!/usr/bin/awk -f
# Insert @PG lines from a saved CRAM header into a SAM/BAM stream.
#
# Usage: awk -v pgfile=<cram_pg.tmp> -f insert_cram_pg_header
#
# For every record stream (SAM/BAM piped via stdin), on the first @PG line:
#   1. Read the saved CRAM @PG lines (pgfile), printing them unchanged.
#      Build a set of their IDs so we can detect name collisions.
#   2. If the current aligner's @PG ID collides with one in the CRAM history,
#      rename the aligner entry by appending .1 (.2, .3, …) until unique –
#      matching the samtools deduplication convention.
#   3. Inject a PP tag on the (possibly renamed) aligner @PG, linking it to
#      the last ID from the inserted CRAM chain.
# All subsequent lines are passed through unchanged.
function get_id(record) {
    n = split(record, t, "\t")
    for (i = 1; i <= n; i++)
        if (t[i] ~ /^ID:/)
            return substr(t[i], 4)
    return ""
}

/^@PG/ && !pg_i {
    aid = get_id($0)
    pg_count = 0
    last_id = ""
    while ((getline line < pgfile) > 0) {
        print line
        pg_lines[++pg_count] = line
        if ((id = get_id(line)) != "") {
            seen[id] = 1
            last_id = id
        }
    }
    close(pgfile)
    pg_i = 1
    if (aid in seen) {
        sfx = 1
        while ((aid "." sfx) in seen) sfx++
        new_aid = aid "." sfx
        n = split($0, af, "\t"); out = ""
        for (j=1; j<=n; j++) { f = af[j]; if (f ~ /^ID:/) f = "ID:" new_aid; out = out (j > 1 ? "\t" : "") f }
        $0 = out
    }
    if (last_id != "") { print $0 "\tPP:" last_id; next }
}
{ print }
