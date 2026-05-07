#!/bin/bash
# =============================================================================
# SUPPA Results Comparison Script
# Author: Gricey
# Description: Compares significant splicing events across different SUPPA runs
#              Run this from inside the data/ folder
# Usage: bash compare_runs.sh
# =============================================================================

# =============================================================================
# RUNS TO COMPARE — add or remove folders here
# =============================================================================

RUNS=("output_strict" "output_variable_1nt" "output_variable_3nt" "output_variable_5nt" "output_variable_10nt" "output_variable_20nt")
LABELS=("Strict" "Variable 1nt" "Variable 3nt" "Variable 5nt" "Variable 10nt" "Variable 20nt")

EVENT_TYPES=("SE" "A3" "A5" "MX" "RI" "AF" "AL")

P_THRESHOLD=0.05
DPSI_FILTER=0.1

# =============================================================================
# HELPER FUNCTION — count significant events
# =============================================================================

count_sig() {
    local file=$1
    if [ -f "$file" ]; then
        awk -v p="$P_THRESHOLD" 'NR>1 && $3!="nan" && $3+0 < p && ($2 < 0 ? -$2 : $2) >= 0.1' "$file" | wc -l | tr -d ' '
    else
        echo "N/A"
    fi
}

# =============================================================================
# SECTION 1 — events detected by generateEvents
# =============================================================================

echo ""
echo "============================================================"
echo "SECTION 1 — Events detected by generateEvents (IOE lines)"
echo "============================================================"
echo ""

printf "%-6s" "Event"
for label in "${LABELS[@]}"; do
    printf "%-20s" "$label"
done
echo ""

printf "%-6s" "-----"
for label in "${LABELS[@]}"; do
    printf "%-20s" "--------------------"
done
echo ""

for event in "${EVENT_TYPES[@]}"; do
    printf "%-6s" "$event"
    for run in "${RUNS[@]}"; do
	        ioe_strict="${run}/events/events_${event}_strict.ioe"
        
	ioe_var=$(ls ${run}/events/events_${event}_variable_*.ioe 2>/dev/null | head -1)

        if [ -f "$ioe_strict" ]; then
            count=$(wc -l < "$ioe_strict")
        elif [ -n "$ioe_var" ]; then
            count=$(wc -l < "$ioe_var")
        else
            count="N/A"
        fi
        printf "%-20s" "$count"
    done
    echo ""
done

# =============================================================================
# SECTION 2 — significant events CT8 vs CT20
# =============================================================================

echo ""
echo "============================================================"
echo "SECTION 2 — Significant events CT8 vs CT20 (p < $P_THRESHOLD, |dPSI| >= $DPSI_FILTER)"
echo "============================================================"
echo ""

printf "%-6s" "Event"
for label in "${LABELS[@]}"; do
    printf "%-20s" "$label"
done
echo ""

printf "%-6s" "-----"
for label in "${LABELS[@]}"; do
    printf "%-20s" "--------------------"
done
echo ""

ct8_totals=()
for run in "${RUNS[@]}"; do
    ct8_totals+=("0")
done

for event in "${EVENT_TYPES[@]}"; do
    printf "%-6s" "$event"
    idx=0
    for run in "${RUNS[@]}"; do
        file="${run}/diff/diff_${event}.dpsi.temp.0"
        count=$(count_sig "$file")
        printf "%-20s" "$count"
        if [ "$count" != "N/A" ]; then
            ct8_totals[$idx]=$((${ct8_totals[$idx]} + $count))
        fi
        idx=$((idx + 1))
    done
    echo ""
done

printf "%-6s" "TOTAL"
for total in "${ct8_totals[@]}"; do
    printf "%-20s" "$total"
done
echo ""

# =============================================================================
# SECTION 3 — significant events CT20 vs PerKO
# =============================================================================

echo ""
echo "============================================================"
echo "SECTION 3 — Significant events CT20 vs PerKO (p < $P_THRESHOLD, |dPSI| >= $DPSI_FILTER)"
echo "============================================================"
echo ""

printf "%-6s" "Event"
for label in "${LABELS[@]}"; do
    printf "%-20s" "$label"
done
echo ""

printf "%-6s" "-----"
for label in "${LABELS[@]}"; do
    printf "%-20s" "--------------------"
done
echo ""

perko_totals=()
for run in "${RUNS[@]}"; do
    perko_totals+=("0")
done

for event in "${EVENT_TYPES[@]}"; do
    printf "%-6s" "$event"
    idx=0
    for run in "${RUNS[@]}"; do
        file="${run}/diff/diff_${event}.dpsi.temp.1"
        count=$(count_sig "$file")
        printf "%-20s" "$count"
        if [ "$count" != "N/A" ]; then
            perko_totals[$idx]=$((${perko_totals[$idx]} + $count))
        fi
        idx=$((idx + 1))
    done
    echo ""
done

printf "%-6s" "TOTAL"
for total in "${perko_totals[@]}"; do
    printf "%-20s" "$total"
done
echo ""

# =============================================================================
# SECTION 4 — top 10 most significant SE events per run
# =============================================================================

echo ""
echo "============================================================"
echo "SECTION 4 — Top 10 most significant SE events (CT8 vs CT20)"
echo "============================================================"

for i in "${!RUNS[@]}"; do
    run="${RUNS[$i]}"
    label="${LABELS[$i]}"
    file="${run}/diff/diff_SE.dpsi.temp.0"
    echo ""
    echo "--- ${label} ---"
    if [ -f "$file" ]; then
        printf "%-60s %-15s %-10s\n" "Event" "dPSI" "p-value"
        awk 'NR>1 && $3!="nan" && $3+0 < 0.05 {print $0}' "$file" \
            | sort -k3 -n \
            | head -10 \
            | awk '{printf "%-60s %-15s %-10s\n", $1, $2, $3}'
    else
        echo "  File not found: $file"
    fi
done

echo ""
echo "============================================================"
echo "Analysis complete!"
echo "============================================================"
echo ""
