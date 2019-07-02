#!/bin/bash

# run cactus
# seqfile could come from ../hgsvc/cactus-input.sh
# then put through ../cactus/rmask.sh
# need to source /ebs1/cactus/cactus_env/bin/activate first
#

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <SEQFILE> <HALFILE>\n"
	 printf "Arguments:\n"
	 printf "   SEQFILE: Cactus seqfile containing seqfile\n"
	 printf "   HALFILE: \n"
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
HALFILE="${1}"
shift

if [ -f "${HALFILE}" ]
then
    echo "hal exists, not re-rerunning cactus"
else
    /usr/bin/time -v cactus jobstore_${HALFILE} ${SEQFILE} ${HALFILE}
fi
