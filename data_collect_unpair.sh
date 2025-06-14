#!/bin/bash
set -x
runs=("SRR12841967" "SRR12841966" "SRR12841965")
output_dir="Aneurysm"

mkdir -p ${output_dir}
rm -rf tmp
for run in "${runs[@]}"; do
    mkdir -p tmp
    run_name="${run// /_}"
    touch tmp/${run_name}.fastq
    for acc in $run; do
        fasterq-dump -p -O tmp ${acc}
        if [ "tmp/${acc}.fastq" != "tmp/${run_name}.fastq" ]; then
            cat tmp/${acc}.fastq >> tmp/${run_name}.fastq
        fi
    done
    round_million=$(awk -v lines=$(wc -l < tmp/${run_name}.fastq) \
                      -v bytes=$(stat -c %s tmp/${run_name}.fastq) \
                      'BEGIN {
                          val = (lines / bytes) * (500 * 1024 * 1024) * 5;
                          rounded = int((val + 500000) / 1000000) * 1000000;
                          print rounded
                      }')
    head -n ${round_million} tmp/${run_name}.fastq > tmp/${run_name}_half.fastq
    mv tmp/${run_name}_half.fastq tmp/${run_name}.fastq
    fastp -i tmp/${run_name}.fastq -o tmp/${run_name}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${run_name}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${run_name}.csv
    rm -rf tmp
done