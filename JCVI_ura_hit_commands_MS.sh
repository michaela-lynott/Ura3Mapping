#!/bin/sh
# -*- mode: sh; -*-
# vim: set filetype=sh :
#written by Alex Richter and edited by Michaela Lynott for mapping URA3 insertion points along K.max genome

set -u # no uninitialized
set -e # exit on error ( use "|| true" or wrap in "set +e" "set -e" to block )
set -o pipefail # also exit if error in piped command -- must disable for acceptable errors
IFS=$'\n\t' # spaces don't split items
bindir=$(cd `dirname $0` && pwd)

usage () {
    cat <<_USAGE
This program does ...

_USAGE
}

onexit () {
    : # necessary cleanup steps
}
trap onexit EXIT


# option parsing
while getopts "h" opt; do
    case $opt in
        h)
            usage >&2
            exit 0
        ;;
    esac
done
shift $((OPTIND-1))


ds=$1
fq1=$2
ura3_bam1=ura3_vs_${ds}_1.bam
bwa mem ../reference/ura3.insert.fasta $fq1 | samtools view -F 0x04 -b > $ura3_bam1

flank_fa=${ds}_ura3_5p_flank.fna
$bindir/extract_flanking_from_sam.py $ura3_bam1 > $flank_fa

bam_kmax_map=kmax_vs_flank_${ds}.bam
bwa mem ../reference/combined.fasta $flank_fa | samtools view -F 0x04 -b > $bam_kmax_map

bed_kmax_map=kmax_insertions_${ds}.bed
$bindir/sam_hit_boundaries.py $bam_kmax_map > $bed_kmax_map
