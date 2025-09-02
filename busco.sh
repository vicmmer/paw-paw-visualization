#!/bin/bash

# Set lineage and mode
LINEAGE="embryophyta_odb10"
MODE="protein"
THREADS=30

# Create main results folder
mkdir -p busco_results_new

# Loop over each annotated FASTA file
for fasta in *_annotated.fasta; do
    # Get the base name (e.g., Annona_cherimola)
    BASENAME=$(basename "$fasta" _annotated.fasta)

    # Create subfolder for this run
    OUTDIR="busco_results_new/${BASENAME}"
    mkdir -p "$OUTDIR"

    # Run BUSCO
    echo "Running BUSCO on $fasta..."
    busco -i "$fasta" \
          -o "$BASENAME" \
          -l "$LINEAGE" \
          -m "$MODE" \
          -c "$THREADS" \
          -f
          --out_path "$OUTDIR"

done

echo "All BUSCO analyses complete."
