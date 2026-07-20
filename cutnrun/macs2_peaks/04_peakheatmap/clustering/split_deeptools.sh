#!/usr/bin/env bash
#
# split_bed_by_group.sh
#
# Split a deepTools "sorted regions" BED file into one BED file per cluster,
# keyed on the trailing `deepTools_group` column.
#
# Background:
#   plotHeatmap / plotProfile with --outFileSortedRegions emit a BED file that
#   carries a single '#'-prefixed header line and appends a `deepTools_group`
#   field (cluster_1, cluster_2, ...) as the LAST column. This script groups the
#   data rows by that field and writes each group to its own file, copying the
#   header into every output so the split files stay self-describing.
#
# Usage:
#   ./split_bed_by_group.sh input.bed [output_dir] [--strip-group]
#
# Arguments:
#   input.bed       deepTools sorted-regions BED (header line beginning with '#')
#   output_dir      directory to write split files into (default: current dir)
#   --strip-group   drop the trailing deepTools_group column, yielding plain
#                   BED12 that downstream tools (bedtools, HOMER, ...) won't
#                   trip over. Omit to preserve the column verbatim.
#
# Output:
#   <output_dir>/<group>.bed   e.g. cluster_1.bed, cluster_2.bed, ...

set -euo pipefail

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; }

# ---- Parse arguments -------------------------------------------------------
strip_group=0
positional=()
for arg in "$@"; do
    case "$arg" in
        --strip-group) strip_group=1 ;;
        -h|--help)     usage; exit 0 ;;
        -*)            echo "Unknown option: $arg" >&2; exit 1 ;;
        *)             positional+=("$arg") ;;
    esac
done

input="${positional[0]:-}"
outdir="${positional[1]:-.}"

[[ -n "$input" ]]     || { echo "Error: no input BED given." >&2; usage; exit 1; }
[[ -f "$input" ]]     || { echo "Error: input not found: $input" >&2; exit 1; }
mkdir -p "$outdir"

# ---- Split -----------------------------------------------------------------
# All the real work is in awk: locate the group column by NAME from the header
# (so we never hard-code its position), then stream rows to per-group files.
awk -v outdir="$outdir" -v strip="$strip_group" '
function joined(skip,    i, s) {
    # Re-join the current record with tabs, omitting field number `skip`.
    s = ""
    for (i = 1; i <= NF; i++) {
        if (i == skip) continue
        s = (s == "" ? $i : s OFS $i)
    }
    return s
}
function sanitize(g) {
    # Keep filenames safe: collapse anything unusual to underscore.
    gsub(/[^A-Za-z0-9._-]/, "_", g)
    return g
}
BEGIN { FS = OFS = "\t" }

# --- Header line (deepTools writes exactly one, prefixed with #) ------------
NR == 1 && /^#/ {
    for (i = 1; i <= NF; i++)
        if ($i == "deepTools_group") groupCol = i
    if (groupCol == 0) groupCol = NF          # fall back to last column
    headerFull  = $0
    headerStrip = joined(groupCol)
    hasHeader   = 1
    next
}

# --- Data lines ------------------------------------------------------------
{
    if (groupCol == 0) groupCol = NF          # also handles a header-less file

    group = $groupCol
    file  = outdir "/" sanitize(group) ".bed"

    # First time we see a group: (re)create its file and write the header.
    if (!(group in seen)) {
        seen[group] = 1
        n++
        if (hasHeader) {
            hdr = (strip == "1") ? headerStrip : headerFull
            print hdr > file                  # ">" truncates on first open,
        }                                     #     appends thereafter
    }

    row = (strip == "1") ? joined(groupCol) : $0
    print row > file
    counts[group]++
}

END {
    printf("Split into %d file(s) in %s/:\n", n, outdir) > "/dev/stderr"
    for (g in counts)
        printf("  %-20s %6d regions\n", g ".bed", counts[g]) > "/dev/stderr"
}
' "$input"