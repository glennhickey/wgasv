#!/bin/bash
#
# make an input file for cactus given some fastas
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

if [[ "$#" -lt "1" ]]; then
    # Too few arguments
    usage
fi

CHROM="${1}"
shift

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz
gzip -d panTro6.fa.gz
samtools faidx panTro6.fa ${CHROM} > panTro6_${CHROM}.fa
rm panTro6.fa

echo "(((((HG00514.0:0.001,HG00514.1:0.001):0.002,(HG00733.0:0.001,HG00733.1:0.001):0.002):0.002,(NA19240.00:0.001,NA19240.1:0.001):0.004):0.004,hg38:0.001):0.2,chimp:0.2);"
echo "hg38 $(pwd)/hg38_${CHROM}.fa"
echo "HG00514.0 $(pwd)/HG00514_${CHROM}_0.fa"
echo "HG00514.1 $(pwd)/HG00514_${CHROM}_1.fa"
echo "HG00733.0 $(pwd)/HG00733_${CHROM}_0.fa"
echo "HG00733.1 $(pwd)/HG00733_${CHROM}_1.fa"
echo "NA19240.0 $(pwd)/NA19240_${CHROM}_0.fa"
echo "NA19240.1 $(pwd)/NA19240_${CHROM}_1.fa"
echo "chimp $(pwd)/panTro6_${CHROM}.fa"

