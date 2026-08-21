#!/bin/bash
# =============================================================================
# igv_filter.sh
# Author: Gricey
#
# Picks IGV candidate genes by gene-level event *count*, not by dPSI/p-value
# (that's what get_igv_candidates.sh does). The idea: a gene with very few
# events in strict mode but several more in variable mode is a clean visual
# case for IGV — few enough events that the screenshot isn't cluttered, but
# a clear jump so the boundary-relaxation effect is obvious to eyeball.
#
# For each gene, compares:
#   strict = number of events in output_strict
#   var1   = number of events in output_variable_1nt
#   diff   = var1 - strict (how many events were "gained" by relaxing 1nt)
#
# Only keeps genes where strict count <= STRICT_MAX (still simple in strict
# mode) AND diff is within [DIFF_MIN, DIFF_MAX] (a modest, inspectable gain).
#
# USAGE:
#   bash igv_filter.sh <EVENT> <STRICT_MAX> <DIFF_MIN> <DIFF_MAX> [TOP_N]
#
# Examples:
#   bash igv_filter.sh RI 3 1 4        # RI events, strict<=3 genes, gained 1-4
#   bash igv_filter.sh A3 2 1 2 30     # A3 events, top 30 instead of default 20
#   bash igv_filter.sh SE 1 1 1        # SE events, strict has exactly 1, gained exactly 1
# ============================================================

EVENT=${1}
STRICT_MAX=${2}
DIFF_MIN=${3}
DIFF_MAX=${4}
TOP_N=${5:-20}        # optional, defaults to 20

# --- validation ---
if [ -z "$EVENT" ] || [ -z "$STRICT_MAX" ] || [ -z "$DIFF_MIN" ] || [ -z "$DIFF_MAX" ]; then
  echo "Usage: bash igv_filter.sh <EVENT> <STRICT_MAX> <DIFF_MIN> <DIFF_MAX> [TOP_N]"
  echo "Example: bash igv_filter.sh RI 3 1 4"
  exit 1
fi

# Resolve data/ relative to this script so it works from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/../data"

STRICT_IOE="${BASE_DIR}/output_strict/events/events_${EVENT}_strict.ioe"
VAR_IOE="${BASE_DIR}/output_variable_1nt/events/events_${EVENT}_variable_1.ioe"

if [ ! -f "$STRICT_IOE" ] || [ ! -f "$VAR_IOE" ]; then
  echo "ERROR: IOE file(s) not found for event type '$EVENT'."
  echo "  Expected: $STRICT_IOE"
  echo "  Expected: $VAR_IOE"
  exit 1
fi

# Count events per gene (col 2 of the IOE file = gene_id) in each mode.
awk 'NR>1{print $2}' "$STRICT_IOE" | sort | uniq -c | awk '{print $2, $1}' > /tmp/s_${EVENT}.txt
awk 'NR>1{print $2}' "$VAR_IOE"    | sort | uniq -c | awk '{print $2, $1}' > /tmp/v_${EVENT}.txt

echo ""
echo "  EVENT: ${EVENT} | strict_max: ${STRICT_MAX} | diff: ${DIFF_MIN}-${DIFF_MAX} | top: ${TOP_N}"
echo ""
echo "  diff  gene_id                  strict  var1"
echo "  ----  -----------------------  ------  ----"

# Join on gene_id (present in both files), compute diff = var1 - strict,
# apply the strict_max / diff_min / diff_max filters, then take the top N
# sorted by diff (ascending) then variable count.
join /tmp/s_${EVENT}.txt /tmp/v_${EVENT}.txt \
  | awk -v smax=$STRICT_MAX -v dmin=$DIFF_MIN -v dmax=$DIFF_MAX \
    '{
      diff = $3 - $2
      if ($2 <= smax && diff >= dmin && diff <= dmax)
        printf "  %-4d  %-23s  %-6d  %d\n", diff, $1, $2, $3
    }' \
  | sort -k1,1n -k3,3n \
  | head -$TOP_N

echo ""
