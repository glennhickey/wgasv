#!/bin/bash
#
# fix up a vg deconstruct output so we can compare it using toil-vg vcfeval
#                                                                                                                                                                                    

usage() {
    # Print usage to stderr                                                                                                                                                          
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <VCF> <SAMPLE> <CHROM>\n"
    printf "Arguments:\n"
	 printf "   VCF: Input VCF\n"
	 printf "   SAMPLE:\n"
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

if [[ "$#" -lt "3" ]]; then
    # Too few arguments
    usage
fi

VCF="${1}"
shift
SAMPLE="${1}"
shift
CHROM="${1}"
shift

# extract the two haplotypes for the sample
# (note: can't get GT filter (ex  --exclude 'GT="0" || GT="."') to work below, so using grep)
vcfsort ${VCF} | bcftools view - -as ${SAMPLE}_${CHROM}_0 | grep -Fv 0:${SAMPLE} | grep -Fv .:. | bcftools annotate - -x INFO/AC,INFO/AN  > ${VCF}.0.vcf
vcfsort ${VCF} | bcftools view - -as ${SAMPLE}_${CHROM}_1 | grep -Fv 0:${SAMPLE} | grep -Fv .:. | bcftools annotate - -x INFO/AC,INFO/AN  > ${VCF}.1.vcf

# merge them into one
# (this script is in: wgasv/hgsvc/sv-genotyping-paper/human/hgsvc)
hap-merge.py  ${VCF}.0.vcf ${VCF}.1.vcf
