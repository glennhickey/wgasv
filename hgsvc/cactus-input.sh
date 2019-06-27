#!/bin/bash
#
# make an input file for cactus given some fastas
#                                                                                                                                                                                    

usage() {
    # Print usage to stderr                                                                                                                                                          
    exec 1>&2
    printf "Usage: $0 [OPTIONS] \n"
         printf "Arguments:\n"
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

if [[ "$#" -lt "0" ]]; then
    # Too few arguments
    usage
fi


echo "((((HG00514.0:0.01,HG00514.1:0.01):0.02,(HG00733.0:0.01,HG00733.1:0.01):0.02):0.02,(NA19240.0:0.01,NA19240.1:0.01):0.04):0.04,hg38:0.01)"
echo "hg38 $(pwd)/hg38_chr21.fa"
echo "HG00514.0 $(pwd)/HG00514_chr21_0.fa"
echo "HG00514.1 $(pwd)/HG00514_chr21_1.fa"
echo "HG00733.0 $(pwd)/HG00533_chr21_0.fa"
echo "HG00733.1 $(pwd)/HG00533_chr21_1.fa"
echo "NA19240.0 $(pwd)/NA19240_chr21_0.fa"
echo "NA19240.1 $(pwd)/NA19240_chr21_1.fa"

