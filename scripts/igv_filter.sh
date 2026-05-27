#!/bin/bash

# ============================================================
# IGV CANDIDATE FILTER — argument-based version
# Usage:
#   bash igv_filter.sh <EVENT> <STRICT_MAX> <DIFF_MIN> <DIFF_MAX> [TOP_N]
#
# Examples:
#   bash igv_filter.sh RI 3 1 4
#   bash igv_filter.sh A3 2 1 2 30
#   bash igv_filter.sh SE 1 1 1
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

STRICT_IOE="/Users/gricey/Desktop/Internship/data/output_strict/events/events_${EVENT}_strict.ioe"
VAR_IOE="/Users/gricey/Desktop/Internship/data/output_variable_1nt/events/events_${EVENT}_variable_1.ioe"

awk 'NR>1{print $2}' "$STRICT_IOE" | sort | uniq -c | awk '{print $2, $1}' > /tmp/s_${EVENT}.txt
awk 'NR>1{print $2}' "$VAR_IOE"    | sort | uniq -c | awk '{print $2, $1}' > /tmp/v_${EVENT}.txt

echo ""
echo "  EVENT: ${EVENT} | strict_max: ${STRICT_MAX} | diff: ${DIFF_MIN}-${DIFF_MAX} | top: ${TOP_N}"
echo ""
echo "  diff  gene_id                  strict  var1"
echo "  ----  -----------------------  ------  ----"

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
