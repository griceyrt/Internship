#!/usr/bin/env bash
# =============================================================================
# get_igv_candidates.sh
# Filters SUPPA diffSplice output and returns IGV-ready coordinates
#
# USAGE:
#   bash scripts/get_igv_candidates.sh
# =============================================================================

# ── CHANGE THESE PARAMETERS ──────────────────────────────────────────────────
RUN="output_strict"
EVENT_TYPE="SE"
COMPARISON="CT8"
DPSI_THRESH=0.1
PVAL_THRESH=0.05
TOP_N=10
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/../data"

GTF="${BASE_DIR}/transcriptome_productivity.gtf"
DIFF_DIR="${BASE_DIR}/${RUN}/diff"

# Pick temp.0 (CT8 vs CT20) or temp.1 (CT20 vs PerKO)
if [ "$COMPARISON" == "CT8" ]; then
    TEMP_FILE="${DIFF_DIR}/diff_${EVENT_TYPE}.dpsi.temp.0"
else
    TEMP_FILE="${DIFF_DIR}/diff_${EVENT_TYPE}.dpsi.temp.1"
fi

# Special case: SE strict CT8 vs CT20 has a pre-cleaned file
if [ "$EVENT_TYPE" == "SE" ] && [ "$COMPARISON" == "CT8" ] && [ "$RUN" == "output_strict" ]; then
    CLEAN="${DIFF_DIR}/diff_SE_clean.dpsi"
    if [ -f "$CLEAN" ]; then
        TEMP_FILE="$CLEAN"
    fi
fi

echo "============================================================"
echo "  SUPPA IGV Candidate Finder"
echo "  Run:        $RUN"
echo "  Event:      $EVENT_TYPE"
echo "  Comparison: $COMPARISON vs $([ "$COMPARISON" == "CT8" ] && echo "CT20" || echo "PerKO")"
echo "  |dPSI| >=  $DPSI_THRESH"
echo "  p-val <=   $PVAL_THRESH"
echo "  File:       $TEMP_FILE"
echo "============================================================"
echo ""

if [ ! -f "$TEMP_FILE" ]; then
    echo "ERROR: File not found: $TEMP_FILE"
    echo "Check that RUN and EVENT_TYPE are correct."
    exit 1
fi

# Filter and sort
RESULTS=$(awk -v dpsi="$DPSI_THRESH" -v pval="$PVAL_THRESH" '
    NR>1 && $3!="nan" && $3+0 < pval && ($2 < 0 ? -$2 : $2) >= dpsi {print $0}
' "$TEMP_FILE" | sort -k3 -g | head -n "$TOP_N")

if [ -z "$RESULTS" ]; then
    echo "No significant events found with these thresholds."
    exit 0
fi

echo "── Top $TOP_N significant $EVENT_TYPE events ($RUN) ─────────────────────"
echo ""

# Parse coordinates per event type and print IGV-ready info
echo "$RESULTS" | while IFS=$'\t' read -r event_id dpsi pval; do

    gene_id=$(echo "$event_id" | cut -d';' -f1)
    chr=$(echo "$event_id" | grep -oE ':[0-9XY]+:' | head -1 | tr -d ':')

    # Get gene name from GTF
    gene_name=$(grep -w "$gene_id" "$GTF" 2>/dev/null \
        | awk '$3=="gene" {
            for(i=1;i<=NF;i++) {
                if($i=="gene_name") {gsub(/"/, "", $(i+1)); print $(i+1)}
            }
        }' | head -1)

    # Parse red region based on event type
    case "$EVENT_TYPE" in
        SE)
            coordA2=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | head -1 | cut -d'-' -f2)
            coordB1=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | tail -1 | cut -d'-' -f1)
            red_start=$coordA2
            red_end=$coordB1
            ;;
        RI)
            retained=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | head -1)
            red_start=$(echo "$retained" | cut -d'-' -f1)
            red_end=$(echo "$retained"   | cut -d'-' -f2)
            ;;
        A3)
            s2=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | head -1 | cut -d'-' -f1)
            s3=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | tail -1 | cut -d'-' -f1)
            if [ "$s3" -lt "$s2" ]; then
                red_start=$s3; red_end=$s2
            else
                red_start=$s2; red_end=$s3
            fi
            ;;
        A5)
            e2=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | head -1 | cut -d'-' -f2)
            e1=$(echo "$event_id" | grep -oE '[0-9]+-[0-9]+' | tail -1 | cut -d'-' -f2)
            if [ "$e1" -lt "$e2" ]; then
                red_start=$e1; red_end=$e2
            else
                red_start=$e2; red_end=$e1
            fi
            ;;
        *)
            red_start="?"
            red_end="?"
            ;;
    esac

    # Calculate size and zoom window
    if [ "$red_start" != "?" ] && [ "$red_end" != "?" ]; then
        size=$((red_end - red_start))
        zoom_start=$((red_start - 500))
        zoom_end=$((red_end + 500))
        zoom="chr${chr}:${zoom_start}-${zoom_end}"
    else
        size="?"
        zoom="n/a"
    fi

    echo "  Gene:        ${gene_name:-[not found in GTF]}  ($gene_id)"
    echo "  dPSI:        $dpsi"
    echo "  p-value:     $pval"
    echo "  IGV search:  ${gene_name:-$gene_id}"
    echo "  Red region:  chr${chr}:${red_start}-${red_end}  (${size} bp)"
    echo "  Zoom to:     $zoom"
    echo "  ──────────────────────────────────────────────────────"
done
