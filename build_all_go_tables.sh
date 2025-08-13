#File to filter interpro scan output into the right format for topGO script from plantbeast https://github.com/PlantBeast/genefamilies/blob/main/TopGO.R 
#!/bin/bash
# Run in the directory that has output tsvs from intepro like OGxxxxx_cleaned.tsv

GO_COL=14
OUT_EXT="all_go_extended.tsv"   # full merge: all rows (no comments), OG appended, transpos removed
OUT_GO="all_go.tsv"             # only OG + GO columns (rows with actual GO terms)
OUT_GO_UNIQ="all_go_unique.tsv" # deduped OG+GO pairs

# start fresh
: > "$OUT_EXT"
: > "$OUT_GO"

# 1–3) Add OG column per file, concatenate all, drop any row containing "transpos" (case-insensitive)
for f in *_cleaned.tsv; do
  [ -f "$f" ] || continue
  og="${f%%_cleaned.tsv}"                # part before _cleaned.tsv
  awk -v og="$og" 'BEGIN{FS=OFS="\t"}
    $0 !~ /^#/ && tolower($0) !~ /transpos/ { print $0, og }' "$f" >> "$OUT_EXT"
done

# 4) Keep only: OG (last column) + GO column (GO_COL), and only if GO column has real GO terms
awk -v c="$GO_COL" 'BEGIN{FS=OFS="\t"}
  $c != "-" && index($c,"GO:") { print $NF, $c }' "$OUT_EXT" > "$OUT_GO"

# 5) Deduplicate OG+GO pairs
LC_ALL=C sort -u "$OUT_GO" > "$OUT_GO_UNIQ"

# 6) Clean GO column: drop "(InterPro)" (and "(PANTHER)" if present), and change '|' to ';'
#sed -E 's/\((InterPro|PANTHER)\)//g; s/\|/;/g' "$OUT_GO_UNIQ" > "$OUT_GO_UNIQ_CLEAN"

awk 'BEGIN{FS=OFS="\t"}
{
  gsub(/\r/,"",$2);                             # strip any Windows CRs in col2
  gsub(/\(InterPro\)|\(PANTHER\)/,"",$2);       # drop source tags
  gsub(/\|/,";",$2);                            # change separators
  print $1, $2
}' all_go_unique.tsv > all_go_unique_clean.tsv
