#!/bin/bash
set -x
accessions=("SRR24084826" "SRR26788820" "SRR26788831" "SRR25437734" "SRR25437738" "SRR24938754" "SRR24938771" "SRR21583662" "ERR7425421" "ERR7425407" "ERR7425402" "ERR7425389" "ERR7425377" "ERR10030092" "SRR19429966" "SRR15539282" "SRR13606063" "SRR11793803" "SRR11793801" "SRR11124505")
output_dir="Pancreatic_Ductal_Adenocarcinoma"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done