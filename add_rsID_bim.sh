#!/bin/bash

# Usage: ./map_rsids.sh file1.txt input.bim output.bim

file1=$1
bim=$2
out=$3

awk '
BEGIN { OFS="\t" }

# -------- Read first file and build mapping --------
FNR==NR {
    # Data comes in blocks of 4 fields:
    # chr:pos  chr  pos  rsID
    for (i=1; i<=NF; i+=4) {
        split($i, a, ":")   # a[1]=chr, a[2]=pos
        key = a[1] "_" a[2]
        rsid = $(i+3)
        map[key] = rsid
    }
    next
}

# -------- Process BIM file --------
{
    chr = $1
    pos = $4
    key = chr "_" pos

    if (key in map) {
        $2 = map[key]   # replace rsID column
    }

    print
}
' "$file1" "$bim" > "$out"
