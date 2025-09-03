#!/bin/bash

# Input folder containing .fa files
input_dir="/home/vmartinez/september2_pipeline_rerun/orthofinder_output/Results_Sep02/Orthogroup_Sequences"

# Directory where this script is located
script_dir="$(cd "$(dirname "$0")" && pwd)"

# Output folder relative to the script directory
output_dir="${script_dir}/interproscan_output"
mkdir -p "$output_dir"

# Path to InterProScan
ips_path="/home/vmartinez/my_interproscan/interproscan-5.75-106.0/interproscan.sh"

# Loop over all .fa files in input_dir
for file in "$input_dir"/*.fa; do
    # Skip if no .fa files are found
    [ -e "$file" ] || { echo "No .fa files found in $input_dir"; exit 1; }

    basename=$(basename "$file" .fa)
    outfile="${output_dir}/${basename}.tsv"

    echo "Running InterProScan on: $file"

    "$ips_path" \
        -i "$file" \
        -f tsv \
        -appl Pfam,PANTHER \
        --iprlookup \
        --goterms \
        -cpu 1 \
        -o "$outfile"
done

#Note: run inside a screen: 
#screen -L -S interproscan 
#find orthofinder_output/Results_Sep02/Orthogroup_Sequences \
#  -maxdepth 1 -type f -name '*.fa' -print0 \
# | parallel -0 -j 40 --eta --verbose ./run_interproscan.sh {}
