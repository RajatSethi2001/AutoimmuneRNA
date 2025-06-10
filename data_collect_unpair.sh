#!/bin/bash
set -x
runs=("SRR24432560" "SRR24432561" "SRR24432562" "SRR24432563" "SRR24432564" "SRR24432565" "SRR24432566" "SRR24432567" "SRR24432568" "SRR24432569" "SRR24432570" "SRR24432572" "SRR24432573" "SRR24432574" "SRR24432575" "SRR24432576" "SRR24432577" "SRR24432578" "SRR24432579" "SRR24432580" "SRR24432581" "SRR24432582" "SRR24432583" "SRR24432584" "SRR24432586" "SRR24432587" "SRR24432588" "SRR24432589" "SRR24432590" "SRR24432592")
output_dir="Scleroderma"

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
    fastp -i tmp/${run_name}.fastq -o tmp/${run_name}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${run_name}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${run_name}.csv
    rm -rf tmp
done