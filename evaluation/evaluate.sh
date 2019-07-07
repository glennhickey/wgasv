#!/bin/bash
#
# compare a graph back to the vcf
#                                                                                                                                                                                    

usage() {
    # Print usage to stderr                                                                                                                                                          
    exec 1>&2
    printf "Usage: $0 [OPTIONS] <GRAPH> <VCF> <REF> <ALT> <SAMPLE> \n"
    printf "Arguments:\n"
	 printf "   GRAPH: Input graph, ex from cactus of seqwish\n"
	 printf "   VCF:   Baseline to compare back to\n"
	 printf "   REF:   Reference path prefix in GRAPH\n"
	 printf "   ALT:   Alt path prefix(es) in GRAPH\n"
	 printf "   SAMPLE:Sample name in VCF\n" 
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

if [[ "$#" -lt "5" ]]; then
    # Too few arguments
    usage
fi

GRAPH="${1}"
shift
VCF="${1}"
shift
REF="${1}"
shift
ALT="${1}"
shift
SAMLPE="${1}"

printf "Need to run in toil-vg virtualenv with vg's .source_me.sh sourced\n"

set +x

# extract a VCF from our graph
OUT_NAME="$(basename $GRAPH .vg).${SAMPLE}"
OUT_VCF="${OUT_NAME}.vcf"

if [ -f ${OUT_VCF} ]
then
    echo "vcf exists, skipping deconstruct"
else
    vg deconstruct ${GRAPH} -P ${REF} -A ${ALT} -e > ${OUT_VCF}
fi
       
# extract the two haplotypes for the sample
# (note: can't get GT filter (ex  --exclude 'GT="0" || GT="."') to work below, so using grep)
vcfsort ${OUT_VCF} | bcftools view - -as ${SAMPLE}_${CHROM}_0 | grep -Fv 0:${SAMPLE} | grep -Fv .:. | bcftools annotate - -x INFO/AC,INFO/AN  > ${OUT_VCF}.0.vcf
vcfsort ${OUT_VCF} | bcftools view - -as ${SAMPLE}_${CHROM}_1 | grep -Fv 0:${SAMPLE} | grep -Fv .:. | bcftools annotate - -x INFO/AC,INFO/AN  > ${OUT_VCF}.1.vcf

# merge them into one
# (this script is in: wgasv/hgsvc/sv-genotyping-paper/human/hgsvc)
wget -nc https://raw.githubusercontent.com/vgteam/sv-genotyping-paper/master/human/hgsvc/hap-merge.py
chmod u+x ./hap-merge.py
./hap-merge.py  ${OUT_VCF}.0.vcf ${OUT_VCF}.1.vcf | bgzip > ${OUT_VCF}.gz
tabix -f -p vcf ${OUT_VCF}.vcf.gz

# extract the sample from the baseline graph
bcftools view -a -s $SAMPLE $VCF -O z --exclude 'GT="0" || GT="0|0" || GT="." || GT=./. || GT="0/0"' > baseline_${SAMPLE}.vcf.gz
tabix -f -p vcf baseline_${SAMPLE}.vcf.gz

# do the toil-vg comparison
rm -rf js eval-${SAMPLE}
toil-vg vcfeval ./js ./eval-${SAMPLE} --sveval --call_vcf ${OUT_VCF}.gz --vcfeval_baseline ${VCF} --min_sv_len 50 --realTimeLogging

# redo the comparison, clipping to the non-repeat regions
wget -nc https://github.com/vgteam/sv-genotyping-paper/blob/master/human/sveval/hg38_non_repeats.bed.gz?raw=true
gzip -dc  hg38_non_repeats.bed.gz > hg38_non_repeats.bed
toil-vg vcfeval ./js ./eval-clip-${SAMPLE} --sveval --call_vcf ${OUT_VCF}.gz --vcfeval_baseline ${VCF} --min_sv_len 50 --realTimeLogging --vcfeval_bed_regions ./hg38_non_repeats.bed



