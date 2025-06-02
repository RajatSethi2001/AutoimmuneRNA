#!/bin/bash
set -x
accessions=("SRR20354605" "SRR20354594" "SRR20354607" "SRR19994439" "SRR19994442" "SRR19994451" "SRR19720421" "SRR19720426" "ERR6929657" "SRR14084332" "SRR14084335" "SRR13895092" "SRR13895087" "SRR13895080" "SRR10173291" "SRR10173294" "SRR10173274" "SRR9956129" "SRR9956138" "SRR9956125" "SRR7949573" "SRR7949568" "SRR7949575" "SRR7170916" "SRR7170908" "SRR6514306" "SRR6514297" "SRR4429149" "ERR1338017" "ERR1338013" "ERR833258" "ERR833252" "ERR833239" "SRR28920066")
output_dir="Parkinson"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done