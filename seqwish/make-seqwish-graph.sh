#!/bin/bash

#
# seqwish to make a graph from a fasta
#

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] FASTA\n"
	 printf "Arguments:\n"
	 printf "   FASTA: Fasta file with all the sequences\n"
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

if [[ "$#" -lt "1" ]]; then
    # Too few arguments
    usage
fi

set -x

FASTA=$1
shift

NAME="${FASTA.*}"

if [ -f "${NAME}.paf" ]
then
    echo "PAF xists, skipping minimap2"
else
    /usr/bin/time -v minimap2 ${FASTA} ${FASTA} -c -X -t $(getconf _NPROCESSORS_ONLN) -x asm20 > ${NAME}.paf 2> TIME.${NAME}.paf
    rm -f ${NAME}.gfa
fi

if [ -f "${NAME}.gfa" ]
then
    echo "GFA exists, skipping seqwish"
else
    /usr/bin/time -v seqwish -s ${FASTA} -p ${NAME}.paf -b ./swtemp -t $(getconf _NPROCESSORS_ONLN) -g ${NAME}.gfa 2> TIME.${NAME}.gfa
    rm -f ${NAME}.vg
fi

if [ -f "${NAME}.vg" ]
then
    echo "vg exists, skipping view"
else
    /usr/bin/time -v vg view -Fv ${NAME}.gfa | vg mod -X 32 - | vg ids -s - > ${NAME}.seqwish.vg 2> TIME.${NAME}.seqwish.vg
fi

