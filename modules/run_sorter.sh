#!/bin/bash
set -e
INPUT_FASTQ=$1
OUTPUT_DIR=$2
MIN_LEN=$3
MAX_LEN=$4
CPUS=$5

export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export TMPDIR=/tmp

python /opt/amplicon_sorter/amplicon_sorter.py \
    -i "$INPUT_FASTQ" \
    -o "$OUTPUT_DIR" \
    -min "$MIN_LEN" \
    -max "$MAX_LEN" \
    -np "$CPUS"