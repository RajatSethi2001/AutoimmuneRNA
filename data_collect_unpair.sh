#!/bin/bash
set -x
accessions=("SRR18185318" "SRR18185286" "SRR18185273" "SRR25660961" "SRR25660964" "SRR25660969" "SRR25660980" "SRR24968372" "SRR24968377" "SRR24968381" "SRR24968392" "SRR24968420" "SRR18185268" "SRR18185278" "SRR18185313" "SRR18185322" "SRR8928991" "SRR8928984" "SRR8928976" "SRR8928970" "SRR3714753" "SRR3714698" "SRR3714650" "SRR3493783" "SRR3493777" "ERR248700" "SRR11991289" "SRR11991284" "SRR26827140" "SRR26827153" "SRR26827151" "SRR24968395" "SRR24968404" "SRR8928967")
output_dir="Ulcerative_Colitis"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp -e 3 ${accession}
    nice fastp -i tmp/${accession}.fastq -o tmp/${accession}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${accession}_trim.fq -p 3 -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done