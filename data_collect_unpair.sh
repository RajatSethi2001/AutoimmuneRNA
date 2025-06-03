#!/bin/bash
set -x
accessions=("SRR12850851" "SRR12850849" "SRR12850840" "SRR12158159" "SRR12158156" "SRR12158154" "SRR10965005" "SRR10965001" "SRR10964991" "SRR10257347" "SRR10257345")
output_dir="Alzheimer"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done