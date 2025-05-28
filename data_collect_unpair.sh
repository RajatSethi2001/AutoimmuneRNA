#!/bin/bash
set -ex
accessions=("SRR18936667" "SRR18936668" "SRR18936672" "SRR18313358" "SRR18313364" "SRR18313398" "SRR18313413" "SRR17155463" "SRR13822264" "SRR13276279" "SRR13276276" "SRR12539734" "SRR12539730" "SRR11783201" "ERR3289166" "ERR2179082")
output_dir="Rheumatoid_Arthritis"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp -e 3 ${accession}
    fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -p 3 -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done