
#!/bin/bash
set -e
accessions=("SRR30764753")
output_dir="Lupus"

mkdir -p ${output_dir}
for accession in "${accessions[@]}"; do
    fasterq-dump -p -O tmp -e 3 ${accession}
    fastp -i tmp/${accession}_1.fastq -I tmp/${accession}_2.fastq -o tmp/${accession}_trim_1.fq -O tmp/${accession}_trim_2.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -1 tmp/${accession}_trim_1.fq -2 tmp/${accession}_trim_2.fq -p 3 -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${accession}.csv
    rm -rf tmp
done