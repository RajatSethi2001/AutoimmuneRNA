#!/bin/bash
set -x
accessions=("SRR9835660" "SRR9835663" "SRR9835656" "SRR9835787" "ERR1993134" "ERR1993135" "ERR1993127" "SRR5659437" "SRR5659434" "SRR5176867")
output_dir="Crohn"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp -e 3 ${accession}
    nice fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -p 3 -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done