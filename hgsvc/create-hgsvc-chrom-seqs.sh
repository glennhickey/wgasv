#!/bin/bash

#
# uses toil-vg (locally) to construct a HGSVC graph for a given chromosome
# haplotypes then get extracted from the graph via the gbwt into fasta files
#

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <CHROM>\n"
	 printf "Arguments:\n"
	 printf "   CHROM: Name of chromosome. ex chr21\n"
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

if [[ "$#" -lt "1" ]]; then
    # Too few arguments
    usage
fi

set -x

CHROM="${1}"
shift

#
# Get the HGSVC VCF via the sv paper scripts
#
if [ -f "sv-genotyping-paper/human/hgsvc/HGSVC.haps.vcf.gz" ]
then
	 echo "HGSVC vcf exists, not regenerating"
else
	 git clone https://github.com/vgteam/sv-genotyping-paper.git --recursive
	 pushd .
	 cd sv-genotyping-paper/human/hgsvc
	 ./make-vcf.sh
	 popd
fi

#
# Use toil-vg to make a chromosome graph
#
NAME="HGSVC_${CHROM}"
if [ -f "${NAME}.gbwt" ]
then
	 echo "Graph exists, not reconstructing"
else
	 toil clean jobstore
	 toil-vg construct jobstore . --regions ${CHROM} --fasta sv-genotyping-paper/human/hgsvc/hg38.fa.gz --vcf sv-genotyping-paper/human/hgsvc/HGSVC.haps.vcf.gz --pangenome --realTimeLogging --workDir . --container None --flat_alts --xg_index --gbwt_index --snarls_index --out_name $NAME --gbwt_index_cores $(getconf _NPROCESSORS_ONLN)
fi

#
# Use vg paths to extract our fasta sequences
#
vg paths -x ${NAME}.xg -F -Q $CHROM > hg38_${CHROM}.fa
cat hg38_${CHROM}.fa > hgsvc_${CHROM}.fa
for SAMPLE in HG00514 HG005733 NA19240
do
    vg paths -x ${NAME}.xg -g ${NAME}.gbwt -F -Q  _thread_${SAMPLE}_${CHROM}_0 | sed 's/_thread_//g'  > ${SAMPLE}_${CHROM}_0.fa
    cat ${SAMPLE}_${CHROM}_0.fa >> hgsvc_${CHROM}.fa
    vg paths -x ${NAME}.xg -g ${NAME}.gbwt -F -Q  _thread_${SAMPLE}_${CHROM}_1 | sed 's/_thread_//g' > ${SAMPLE}_${CHROM}_1.fa
    cat ${SAMPLE}_${CHROM}_1.fa >> hgsvc_${CHROM}.fa
done
