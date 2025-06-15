#!/bin/bash
ulimit -v $((4000 * 1024))
set -x
runs=("SRR10401138")
output_dir="Lupus"
file_limit=500000000
compress_factor=5

mkdir -p ${output_dir}
rm -rf tmp
sample_spots=100000
true_file_limit=$(($compress_factor * $file_limit))
for acc in "${runs[@]}"; do
    mkdir -p tmp
    fastq-dump -N 1 -X ${sample_spots} -O tmp ${acc}
    sample_bytes=$(du -b tmp/${acc}.fastq | awk '{sum += $1} END {print sum}')
    total_spots=$(($true_file_limit * $sample_spots / $sample_bytes))
    
    fastq-dump -N 1 -X ${total_spots} -O tmp ${acc} &
    DUMP_PID=$!

    set +x
    while kill -0 $DUMP_PID 2>/dev/null; do
        bytes=$(du -b tmp/${acc}.fastq 2>/dev/null | awk '{sum += $1} END {print sum}')
        percent=$((100 * $bytes / $true_file_limit))
        echo "[`date +%T`] Estimated progress: $percent% ($bytes / $true_file_limit Bytes)"
        sleep 10
    done
    set -x

    fastp -i tmp/${acc}.fastq -o tmp/${acc}_trim.fq -j tmp/fastp.json -h tmp/fastp.html
    salmon quant -i salmon_index -l A -r tmp/${acc}_trim.fq -o tmp --validateMappings
    python tx2gene.py -p tmp/quant.sf -o ${output_dir}/${acc}.csv
    rm -rf tmp
done