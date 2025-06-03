#!/bin/bash
set -x
accessions=()
output_dir="Polycystic_Ovary"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done