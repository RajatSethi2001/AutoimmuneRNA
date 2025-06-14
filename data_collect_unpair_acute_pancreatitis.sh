#!/bin/bash
set -x
runs=("SRR17734547" "SRR17734548" "SRR17734549" "SRR17734550" "SRR17734551" "SRR17734552" "SRR17734553" "SRR17734554" "SRR17734555" "SRR17734556" "SRR17734557" "SRR17734558" "SRR17734559" "SRR17734560" "SRR17734561")
output_dir="Acute_Pancreatitis"

mkdir -p ${output_dir}
rm -rf tmp
for run in "${runs[@]}"; do
    mkdir -p tmp
    run_name="${run// /_}"
    touch tmp/${run_name}.fastq
    for acc in $run; do
        fasterq-dump --spot-limit 14000000 -p -O tmp ${acc}
        if [ "tmp/${acc}.fastq" != "tmp/${run_name}.fastq" ]; then
            cat tmp/${acc}.fastq >> tmp/${run_name}.fastq
        fi
    done
    fastp -i tmp/${run_name}.fastq -o tmp/${run_name}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${run_name}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${run_name}.csv
    rm -rf tmp
done