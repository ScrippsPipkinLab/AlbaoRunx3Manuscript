#!/usr/bin/env bash

INPUT=$1
GENOME=${2:-"genome.fa"}

if [ -z "$INPUT" ]; then
  echo "Usage: $(basename $0) <bed_file> [genome.fa]"
  exit 1
fi

OUT_DIR=$(basename "$INPUT" .bed)

run --cpus=11 dsalbao/homer:5.1 findMotifsGenome.pl \
  bed/$INPUT \
  $GENOME \
  homer_output/$OUT_DIR \
  -size given \
  -p 11 -nomotif