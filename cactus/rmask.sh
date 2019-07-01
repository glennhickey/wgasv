#!/bin/bash

# take a cactus input seqfile and repeatmask every sequence and print out a new seqfile
#

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <SEQFILE> <OUT-DIR>\n"
	 printf "Arguments:\n"
	 printf "   SEQFILE: Cactus seqfile containing seqfile to mask\n"
	 printf "   OUT-DIR: Directory for masked sequences\n"
	 printf "Options:\n"
    exit 1
}

while getopts "h" o; do
    case "${o}" in
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

printf "Need to run in toil-vg virtualenv with vg's .source_me.sh sourced\n"

if [[ "$#" -lt "2" ]]; then
    # Too few arguments
    usage
fi

set -x

SEQFILE="${1}"
shift
OUTDIR="${1}"
shift

mkdir ${OUTDIR}
rm -rf ${OUTDIR}.jobstore

repeatMaskerPipeline.py human ${OUTDIR} $(tail -n +2 ${SEQFILE} | grep -v chimp | grep -v hg38 | awk '{ print $2 }' ) ${OUTDIR}.jobstore
head -1 ${SEQFILE}
tail -n +2 ${SEQFILE} | grep -v chimp | grep -v hg38 | awk -v od="${OUTDIR}" '{ print $1 "\t" od basename $2 ".masked" }'
tail -n +2 ${SEQFILE} | grep chimp
tail -n +2 ${SEQFILE} | grep hg3

