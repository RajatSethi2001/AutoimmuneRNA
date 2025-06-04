#!/bin/bash
set -x
accessions=("ERR2105195" "SRR8439274" "SRR8439276" "SRR8439272" "SRR1510156" "SRR1510150" "SRR1510146" "SRR1510138" "SRR1510128" "SRR309133" "SRR309138")
output_dir="Autism"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done