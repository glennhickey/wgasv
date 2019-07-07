#!/bin/bash

# run cactus
# seqfile could come from ../hgsvc/cactus-input.sh
# then put through ../cactus/rmask.sh
# need to source /ebs1/cactus/cactus_env/bin/activate first
#

CONFIG=0

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <SEQFILE> <HALFILE> <REFNAME>\n"
	 printf "Arguments:\n"
	 printf "   SEQFILE: Cactus seqfile containing seqfile\n"
	 printf "   HALFILE: Output\n"
	 printf "   REFNAME: Reference genome\n"
	 printf "Options:\n"
	 printf "  -c FILE:  Config file\n" 
    exit 1
}

while getopts "hc:" o; do
    case "${o}" in
	c)
	    CONFIG=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

printf "Need to run in toil-vg virtualenv with vg's .source_me.sh sourced\n"

if [[ "$#" -lt "3" ]]; then
    # Too few arguments
    usage
fi

set -x

SEQFILE="${1}"
shift
HALFILE="${1}"
shift
REFNAME="${1}"
shift

set -x

OPTS=""
if [ $CONFIG != 0 ]
then
    OPTS="--config ${CONFIG}"
fi

if [ -f "${HALFILE}" ]
then
    echo "hal exists, not re-rerunning cactus"
else
    /usr/bin/time -v cactus jobstore_${HALFILE} ${SEQFILE} ${HALFILE} ${OPTS} 2> TIME.${HALFILE}
fi

VGFILE="$(basename "${HALFILE}" .hal).vg"
if [ -f "${VGFILE}" ]
then 
    echo "vg exits, skipping hal2vg"
else
    /usr/bin/time -v hal2vg ${HALFILE} --refGenome hg38 --noAncestors --targetGenomes HG00514.0,HG00514.1 --onlySequenceNames > ${VGFILE} 2> TIME.${VGFILE}
fi
