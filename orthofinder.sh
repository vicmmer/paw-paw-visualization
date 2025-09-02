#!/bin/bash

# Set input and output folders
INPUT_DIR="protein_sequences"
OUTPUT_DIR="orthofinder_output"

# Set thread counts
THREADS=20

# Run OrthoFinder
orthofinder -f "$INPUT_DIR" -o "$OUTPUT_DIR" -t "$THREADS" -a 10 -S diamond

# Done!
echo "OrthoFinder finished! Results in: $OUTPUT_DIR/OrthoFinder/Results_*"
