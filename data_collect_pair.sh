#!/bin/bash
set -x
accessions=("SRR28410292" "SRR28410302" "SRR28410307" "SRR28410287" "SRR28410279" "SRR28410278" "SRR28410269" "SRR28410191" "SRR28410262" "SRR28410261" "SRR28410259" "SRR28410256" "SRR28410240" "SRR28410230" "SRR28410229" "SRR28410228" "SRR28410207" "SRR28410130" "SRR28410205" "SRR28410204" "SRR28410158" "SRR28410304" "SRR28410305" "SRR28410320")
output_dir="OCD"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}_1.fastq -I tmp/${accession}_2.fastq -o tmp/${accession}_trim_1.fq -O tmp/${accession}_trim_2.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -1 tmp/${accession}_trim_1.fq -2 tmp/${accession}_trim_2.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done