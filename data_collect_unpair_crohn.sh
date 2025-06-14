#!/bin/bash
set -x
runs=("SRR18719951" "SRR18719949" "SRR18719947" "SRR18719938" "SRR18719920" "SRR18719914" "SRR18719913" "SRR18719912" "SRR18719908" "SRR18719907" "SRR18719881" "SRR18719879" "SRR18719876" "SRR18719875" "SRR18720348" "SRR18720347" "SRR18720345" "SRR18720342" "SRR18720303" "SRR18720299" "SRR18720295" "SRR18720293" "SRR18720289" "SRR18720288" "SRR18719595" "SRR18719594" "SRR18719592" "SRR18719589" "SRR18719551" "SRR18719550" "SRR18719549" "SRR18719548" "SRR18719545" "SRR18719827" "SRR18719824" "SRR18719822" "SRR18719820" "SRR18719818" "SRR18719817" "SRR18720070" "SRR18720068" "SRR18719863" "SRR18719861" "SRR18719857" "SRR18719856" "SRR18719854" "SRR18719853" "SRR18719852" "SRR18719851" "SRR18719849" "SRR18719782" "SRR18719781" "SRR18719780" "SRR18719778" "SRR18719776" "SRR18719775" "SRR18719774" "SRR18719773" "SRR18719772" "SRR18719770" "SRR18719472")
output_dir="Crohns_Disease"

mkdir -p ${output_dir}
rm -rf tmp
for run in "${runs[@]}"; do
    mkdir -p tmp
    run_name="${run// /_}"
    touch tmp/${run_name}.fastq
    for acc in $run; do
        fasterq-dump -p -O tmp ${acc}
        if [ "tmp/${acc}.fastq" != "tmp/${run_name}.fastq" ]; then
            cat tmp/${acc}.fastq >> tmp/${run_name}.fastq
        fi
    done
    fastp -i tmp/${run_name}.fastq -o tmp/${run_name}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${run_name}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${run_name}.csv
    rm -rf tmp
done