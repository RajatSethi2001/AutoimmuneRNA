#!/bin/bash
set -x
runs=("SRR22854021" "SRR22854022" "SRR22854023" "SRR22854024" "SRR22854025" "SRR22854026" "SRR22854027" "SRR22854028" "SRR22854029" "SRR22854030" "SRR22854031" "SRR22854032" "SRR22854033" "SRR22854034" "SRR22854035" "SRR22854036" "SRR22854037" "SRR22854038" "SRR22854039" "SRR22854040" "SRR22854041" "SRR22854042" "SRR22854043" "SRR22854044" "SRR22854045" "SRR22854046" "SRR22854047" "SRR22854048" "SRR22854049" "SRR22854050" "SRR22854051" "SRR22854052" "SRR22854053" "SRR22854054" "SRR22854055" "SRR22854056" "SRR22854057" "SRR22854058" "SRR22854059")
output_dir="MRSA_Bacteremia"

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