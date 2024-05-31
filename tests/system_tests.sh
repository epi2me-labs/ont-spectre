#!/bin/bash
set -euo pipefail

declare -a data=( 
    "hg002_chr18_realvar"
    "hg002_chr19_simulations"
)

tabix_script="import pysam; import sys; pysam.tabix_index(sys.argv[1], preset='vcf', force=True)"

for sample in "${data[@]}"
do
    echo "Processing $sample"

    pfx="tests/data/system_tests/$sample/"
    echo "Data path prefix: $pfx"
    
    # spectre cleanup and run
    rm -rf $pfx/output_spectre
    spectre CNVCaller \
            --bin-size 1000 \
            --coverage $pfx/input/mosdepth/ \
            --snv $pfx/input/wf_snp.vcf.gz \
            --sample-id sample \
            --output-dir $pfx/output_spectre \
            --reference $pfx/input/ref.fna.gz \
            --metadata $pfx/input/metadata.mdr \
            --blacklist $pfx/input/blacklist.bed

    # gziping and indexing for truvari input
    vcf_output="$pfx/output_spectre/sample.vcf"
    python -c "$tabix_script" ${vcf_output}

    # spectre cleanup and run
    rm -rf $pfx/truvari
    truvari bench --base $pfx/expected.vcf.gz \
                  --comp $vcf_output.gz \
                  --output $pfx/truvari \
                  --pctseq 0 \
                  --pctovl 0.8 \
                  --refdist 1000000000 \
                  --sizemax 1000000000 \
                  --chunksize 1000000000

    # get and check f1 score
    f1_score=$(grep '"f1"' $pfx/truvari/summary.json | awk -F ': ' '{print $2}' | tr -d ',')

    if [ "$f1_score" == "1.0" ]; then
        echo -e "F1 Score is equal to 1. \nTest passed!\n"
    else
        echo -e "F1 Score is not equal to 1. \nExiting.\n"
        exit 1
    fi

done


