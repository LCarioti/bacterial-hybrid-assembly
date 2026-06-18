#!/bin/bash

#  Bacterial Hybrid Assembly Pipeline
#  Copyright (C) 2026 LCarioti 
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://gnu.org>.

set -euo pipefail

###############################################################################
# HELP
###############################################################################

show_help() {
cat << EOF

Usage:
  bacterial_pipeline.sh <R1> <R2> <nanopore.fastq.gz> <reference.fasta> <threads> [mode: standard|aggressive]

Options for mode (optional):
  standard   : Keep reads > 2000 bp (Default. Best if looking for small plasmids)
  aggressive : Keep reads > 6000 bp, top 90% (Best for closing chromosomes and rRNA)

EOF
exit 0
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && show_help

# Check the 5 mandatory options (6 is optional)
if [ $# -lt 5 ] || [ $# -gt 6 ]; then
    show_help
fi

###############################################################################
# INPUT & CONFIGURATION
###############################################################################

R1=$(realpath "$1")
R2=$(realpath "$2")
NANOPORE=$(realpath "$3")
REFERENCE=$(realpath "$4")
THREADS="$5"
MODE="${6:-standard}" # If 6 is empty, use standard 

WORKDIR=$(pwd)

# Universal Software Discovery (Dynamic lookup from user PATH)
FASTP=$(which fastp 2>/dev/null || true)
MINIMAP2=$(which minimap2 2>/dev/null || true)
BWA=$(which bwa 2>/dev/null || true)
SAMTOOLS=$(which samtools 2>/dev/null || true)
FILTLONG=$(which filtlong 2>/dev/null || true)
FLYE=$(which flye 2>/dev/null || true)
POLCA=$(which polca.sh 2>/dev/null || which polca 2>/dev/null || true)
QUAST=$(which quast.py 2>/dev/null || which quast 2>/dev/null || true)

# Pre-flight Check: Ensure all 8 tools are globally accessible
missing_tools=0
for tool_name in "FASTP" "MINIMAP2" "BWA" "SAMTOOLS" "FILTLONG" "FLYE" "POLCA" "QUAST"; do
    if [ -z "${!tool_name}" ]; then
        echo "Errore: Il tool richiesto '$tool_name' non è stato trovato nel tuo PATH." >&2
        missing_tools=$((missing_tools + 1))
    fi
done

if [ "$missing_tools" -gt 0 ]; then
    echo "Errore globale: Mancano $missing_tools dipendenze. Controlla il tuo ambiente o il file README." >&2
    exit 1
fi

###############################################################################
# DIRECTORIES & LOGGING
###############################################################################

RAW="raw"
QC="qc"
NANO="nanopore"
ASM="assembly"
REF="reference"
LOGS="logs"
MAPPING="mapping"

mkdir -p "$RAW" "$QC" "$NANO" "$ASM" "$REF" "$LOGS" "$MAPPING"

exec > >(tee "$LOGS/pipeline.log")
exec 2>&1

echo "PIPELINE START - MODE: $MODE"
date

###############################################################################
# LINK INPUTS
###############################################################################

R1_B=$(basename "$R1")
R2_B=$(basename "$R2")
NANO_B=$(basename "$NANOPORE")
REF_B=$(basename "$REFERENCE")

ln -sf "$R1" "$RAW/$R1_B"
ln -sf "$R2" "$RAW/$R2_B"
ln -sf "$NANOPORE" "$RAW/$NANO_B"
ln -sf "$REFERENCE" "$REF/$REF_B"

###############################################################################
# FASTP (ILLUMINA QC)
###############################################################################

echo "STEP 1 - FASTP"
"$FASTP" -i "$RAW/$R1_B" -I "$RAW/$R2_B" \
-o "$QC/illumina_R1.qc.fastq.gz" -O "$QC/illumina_R2.qc.fastq.gz" \
--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
--length_required 50 --thread "$THREADS" -h "$QC/fastp.html" -j "$QC/fastp.json"

###############################################################################
# FILTLONG & FLYE 
###############################################################################

if [ "$MODE" == "aggressive" ]; then
    echo "STEP 2 - FILTLONG (Aggressive Mode)"
    "$FILTLONG" --min_length 6000 --keep_percent 90 "$RAW/$NANO_B" | gzip -c > "$NANO/nano.fastq.gz"

    echo "STEP 3 - FLYE"
    "$FLYE" --nano-hq "$NANO/nano.fastq.gz" --out-dir "$ASM" --threads "$THREADS"

elif [ "$MODE" == "standard" ]; then
    echo "STEP 2 - FILTLONG (Standard Mode)"
    "$FILTLONG" --min_length 2000 --keep_percent 95 "$RAW/$NANO_B" | gzip -c > "$NANO/nano.fastq.gz"

    echo "STEP 3 - FLYE"
    "$FLYE" --nano-hq "$NANO/nano.fastq.gz" --out-dir "$ASM" --threads "$THREADS" 
else
    echo "Error: Options '$MODE' not recognized. Use 'standard' or 'aggressive'." >&2
    exit 1
fi

###############################################################################
# POLCA 
###############################################################################

echo "STEP 4 - POLCA"
(
  cd "$ASM/"
  "$POLCA" -a assembly.fasta \
  -r "$WORKDIR/$QC/illumina_R1.qc.fastq.gz $WORKDIR/$QC/illumina_R2.qc.fastq.gz" \
  --threads "$THREADS"
)

###############################################################################
# QUAST
###############################################################################

echo "STEP 5 - QUAST"
"$QUAST" "$ASM/assembly.fasta.PolcaCorrected.fa" -r "$REF/$REF_B" -o "$ASM/quast" --threads "$THREADS"

###############################################################################
# REMAPPING (ONT & Illumina)
###############################################################################

echo "STEP 6 - ONT MAPPING"
# Mapping of Filtlong-filtered reads to the polished genome
"$MINIMAP2" -t "$THREADS" -ax map-ont "$ASM/assembly.fasta.PolcaCorrected.fa" "$NANO/nano.fastq.gz" | \
"$SAMTOOLS" sort -@ "$THREADS" -m 2G -o "$MAPPING/ont.bam"
"$SAMTOOLS" index "$MAPPING/ont.bam"
"$SAMTOOLS" flagstat "$MAPPING/ont.bam" > "$MAPPING/ont_flagstat.txt"
# Mapping of fastp-filtered reads to the polished genome
echo "STEP 7 - ILLUMINA MAPPING"
"$BWA" index "$ASM/assembly.fasta.PolcaCorrected.fa"
"$BWA" mem -t "$THREADS" "$ASM/assembly.fasta.PolcaCorrected.fa" "$QC/illumina_R1.qc.fastq.gz" "$QC/illumina_R2.qc.fastq.gz" | \
"$SAMTOOLS" sort -@ "$THREADS" -m 2G -o "$MAPPING/illumina.bam"
"$SAMTOOLS" index "$MAPPING/illumina.bam"
"$SAMTOOLS" flagstat "$MAPPING/illumina.bam" > "$MAPPING/illumina_flagstat.txt"

echo "PIPELINE COMPLETED"
date

