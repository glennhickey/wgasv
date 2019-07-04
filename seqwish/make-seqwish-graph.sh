#!/bin/bash

# run seqwish
# needs minimap2, fpa and seqish in the path
#
# minimap takes as input a fasta file for each sequence
# seqish takes as input the paf output of minimap along with a fasta file with all
# the sequences inside (with reference sequence first)
#

FPALEN=10000

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <ALL-FASTA> <OUT-NAME> <FASTA1> <FASTA2> ... <FASTAN>\n"
	 printf "Arguments:\n"
	 printf "   ALL-FASTA: file containing all the sequences in desired order\n"
	 printf "   OUT-NAME:  base name for output files\n"	 
	 printf "   FASTA1-N:  same sequences as above, but split into one file per genome\n"
	 printf "Options:\n"
	 printf "  -l N:       FPA length filter [default=${FPALEN}]\n"
    exit 1
}

while getopts "hl:" o; do
    case "${o}" in
	l)
	    FPALEN=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "4" ]]; then
    # Too few arguments
    usage
fi

set -x

ALL_FASTA_FILE="${1}"
shift
OUT_NAME="${1}"
shift

FASTA_LIST=""
for var in "$@"
do
    FASTA_LIST="${FASTA_LIST} $var"
done

# get the minimap script and set the threads
wget -nc https://raw.githubusercontent.com/ekg/yeast-pangenome/master/pan-minimap2
sed -i 's/-t 20/-t $(getconf _NPROCESSORS_ONLN)/g' pan-minimap2
chmod u+x ./pan-minimap2

if [ -f ${OUT_NAME}.paf ]
then
    echo "PAF exists, skipping minimap2"
else
    # run the pairwise alignment
    /usr/bin/time -v ./pan-minimap2 $FASTA_LIST > ${OUT_NAME}.paf 2> TIME.${OUT_NAME}.paf
    rm -f ${OUT_NAME}_fpa${FPALEN}.paf
    rm -f ${OUT_NAME}_fpa${FPALEN}.gfa
fi    

# do the fpa filter (note: fpa was installed via conda)
cat ${OUT_NAME}.paf | fpa drop -l ${FPALEN} > ${OUT_NAME}_fpa${FPALEN}.paf

mkdir -p work # a work directory that's local, files created here will be deleted when seqwish completes

if [ -f ${OUT_NAME}_fpa${FPALEN}.gfa ]
then
    echo "GFA exits, skipping seqwish"
else
    /usr/bin/time -v seqwish -s ${ALL_FASTA_FILE} -p ${OUT_NAME}_fpa${FPALEN}.paf -t $(getconf _NPROCESSORS_ONLN) -b work/x -g ${OUT_NAME}_fpa${FPALEN}.gfa 2> TIME.${OUT_NAME}_fpa${FPALEN}.gfa
    rm -f ${OUT_NAME}_fpa${FPALEN}.vg
fi

if [ -f ${OUT_NAME}_fpa${FPALEN}.vg ]
then
    echo "VG exists, skipping vg view"
else
    /usr/bin/time -v vg view -Fv ${OUT_NAME}_fpa${FPALEN}.gfa > ${OUT_NAME}_fpa${FPALEN}.vg 2> TIME.${OUT_NAME}_fpa${FPALEN}.vg
fi


