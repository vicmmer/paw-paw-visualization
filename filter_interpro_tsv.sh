#!/usr/bin/env bash
set -euo pipefail

# 1) Get all the *_cleaned.tsv files in this folder
find . -type f -name '*_cleaned.tsv' -print0 | sort -z > /tmp/all_cleaned.zlist

# 2) Check which files are completely empty (0 bytes)
find . -type f -name '*_cleaned.tsv' -size 0 -print0 | sort -z > /tmp/empty.zlist

# 3) Check which files have the word "transpos" (covers transposon/transposase, not case sensitive)
xargs -0 grep -I -i -l -Z 'transpos' < /tmp/all_cleaned.zlist > /tmp/transpos.zlist || true

# 4) Make normal newline lists so I can look at them easily later
tr '\0' '\n' < /tmp/all_cleaned.zlist > all_cleaned_tsv.txt
tr '\0' '\n' < /tmp/empty.zlist       > empty_tsv.txt
tr '\0' '\n' < /tmp/transpos.zlist    > transpos_tsv.txt

# 5) Merge empty + transpos lists into one "bad files" list
sort -u empty_tsv.txt transpos_tsv.txt > exclude_tsv.txt

# 6) Figure out which files are okay to keep = all - bad ones
sort -u all_cleaned_tsv.txt > /tmp/all.sorted.txt
sort -u exclude_tsv.txt     > /tmp/exclude.sorted.txt
comm -23 /tmp/all.sorted.txt /tmp/exclude.sorted.txt > include_tsv.txt

# 6b) Clean up the file names in the keep list (remove leading ./ and _cleaned.tsv)
sed -E 's|^\./||; s|_cleaned\.tsv$||' include_tsv.txt > include_orthogroups.txt
# Do the same for exclude list if you want cleaned orthogroup names there too
sed -E 's|^\./||; s|_cleaned\.tsv$||' exclude_tsv.txt > exclude_orthogroups.txt

# 7) Quick counts
total=$(wc -l < all_cleaned_tsv.txt)
excluded=$(wc -l < exclude_tsv.txt)
kept=$(wc -l < include_tsv.txt)

# 8) Print the stats
echo "Total *_cleaned.tsv files: $total"
echo "No good (empty or transpos*): $excluded"
echo "Good ones to keep: $kept"
