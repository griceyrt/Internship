#!/usr/bin/env bash
# =============================================================================
# compare_event_populations.sh
#
# Answers Kiran's three questions about strict vs variable 1nt IOE catalogues:
#
#   Q1 — Which genes go from strict → MORE events? (gained events, both present)
#   Q2 — Which genes go from 0 → events?           (absent in strict, new in var)
#   Q3 — Which events involve novel transcripts?    (UUID-format IDs, not ENSMUST)
#
# USAGE:
#   bash compare_event_populations.sh
#
# Place in data/ (same level as output_strict/ and output_variable_1nt/)
# =============================================================================

# ── PATHS ────────────────────────────────────────────────────────────────────
IOE_STRICT="/Users/gricey/Desktop/Internship/data/output_strict/events"
IOE_VAR1="/Users/gricey/Desktop/Internship/data/output_variable_1nt/events"
OUTPUT_DIR="/Users/gricey/Desktop/Internship/data/ioe_populations"
mkdir -p "$OUTPUT_DIR"

EVENTS=("A3" "A5" "RI")
SUFFIX_STRICT="strict"
SUFFIX_VAR1="variable_1"

# ── HEADER ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  IOE Population Comparison: Strict vs Variable 1nt"
echo "============================================================"
echo ""

# ── PER EVENT TYPE ────────────────────────────────────────────────────────────
for EVENT in "${EVENTS[@]}"; do

    STRICT_IOE="${IOE_STRICT}/events_${EVENT}_${SUFFIX_STRICT}.ioe"
    VAR1_IOE="${IOE_VAR1}/events_${EVENT}_${SUFFIX_VAR1}.ioe"

    if [ ! -f "$STRICT_IOE" ] || [ ! -f "$VAR1_IOE" ]; then
        echo "  [SKIP] $EVENT — IOE file(s) not found"
        echo ""
        continue
    fi

    echo "────────────────────────────────────────────────────────────"
    echo "  EVENT TYPE: $EVENT"
    echo "────────────────────────────────────────────────────────────"

    # ── Q1: Genes with MORE events in variable (strict > 0, variable > strict) ──
    # Count events per gene in each mode, then join and find genes where var > strict

    # Count events per gene in strict (skip header)
    awk 'NR>1 {count[$2]++} END {for (g in count) print g, count[g]}' \
        "$STRICT_IOE" | sort > /tmp/strict_counts_${EVENT}.txt

    # Count events per gene in variable
    awk 'NR>1 {count[$2]++} END {for (g in count) print g, count[g]}' \
        "$VAR1_IOE" | sort > /tmp/var1_counts_${EVENT}.txt

    # Join on gene_id, keep only genes present in BOTH, where var > strict
    join /tmp/strict_counts_${EVENT}.txt /tmp/var1_counts_${EVENT}.txt \
        | awk '$3 > $2 {diff=$3-$2; print $1, $2, $3, diff}' \
        | sort -k4 -rn \
        > "${OUTPUT_DIR}/${EVENT}_Q1_more_events.txt"

    Q1_COUNT=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q1_more_events.txt")

    echo ""
    echo "  Q1 — Genes gaining events (strict→more, both present): $Q1_COUNT genes"
    echo "       Columns: gene_id | strict_count | var1_count | gained"
    echo "       Top 10:"
    echo "       gene_id                    strict  var1   gained"
    echo "       ─────────────────────────────────────────────────"
    head -10 "${OUTPUT_DIR}/${EVENT}_Q1_more_events.txt" \
        | awk '{printf "       %-30s %4s    %4s   +%s\n", $1, $2, $3, $4}'

    # ── Q2: Genes absent in strict, present in variable (0 → events) ──────────
    # Get gene sets
    awk 'NR>1 {print $2}' "$STRICT_IOE" | sort -u > /tmp/strict_genes_${EVENT}.txt
    awk 'NR>1 {print $2}' "$VAR1_IOE"   | sort -u > /tmp/var1_genes_${EVENT}.txt

    # Set difference: in var1 but NOT in strict
    comm -13 /tmp/strict_genes_${EVENT}.txt /tmp/var1_genes_${EVENT}.txt \
        > "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_genes.txt"

    Q2_GENES=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_genes.txt")

    # Also get the IOE rows for these genes (with event counts)
    join "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_genes.txt" \
         /tmp/var1_counts_${EVENT}.txt \
        | sort -k2 -rn \
        > "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_counts.txt"

    echo ""
    echo "  Q2 — Genes absent in strict, new in variable (0→events): $Q2_GENES genes"
    echo "       Top 10 by number of new events:"
    echo "       gene_id                    var1_events"
    echo "       ─────────────────────────────────────"
    head -10 "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_counts.txt" \
        | awk '{printf "       %-30s %s\n", $1, $2}'

    # ── Q3: Events with novel (UUID-format) transcript IDs ────────────────────
    # Novel transcripts appear in the alternative_transcripts column (col 3)
    # They match UUID format: 8hex-4hex-4hex-4hex-12hex (not ENSMUST...)
    # In variable IOE, col3 = alternative_transcripts, col4 = total_transcripts

    # Extract events where any transcript ID in col3 or col4 is UUID-format
    awk 'NR>1 {
        # Check alternative_transcripts (col3) and total_transcripts (col4)
        if ($4 ~ /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/ ||
            $5 ~ /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/) {
            print $0
        }
    }' "$VAR1_IOE" > "${OUTPUT_DIR}/${EVENT}_Q3_novel_transcripts.ioe"

    Q3_EVENTS=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q3_novel_transcripts.ioe")

    # Count unique genes involved
    Q3_GENES=$(awk '{print $2}' "${OUTPUT_DIR}/${EVENT}_Q3_novel_transcripts.ioe" | sort -u | wc -l)

    echo ""
    echo "  Q3 — Events with novel (UUID) transcripts: $Q3_EVENTS events across $Q3_GENES genes"
    echo "       Top 10 genes by novel event count:"
    awk '{print $2}' "${OUTPUT_DIR}/${EVENT}_Q3_novel_transcripts.ioe" \
        | sort | uniq -c | sort -rn | head -10 \
        | awk '{printf "       %-30s %s novel events\n", $2, $1}'

    echo ""

done

# ── SUMMARY TABLE ─────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  SUMMARY"
echo "============================================================"
printf "  %-6s  %10s  %10s  %10s\n" "Event" "Q1 gained" "Q2 new" "Q3 novel"
printf "  %-6s  %10s  %10s  %10s\n" "------" "---------" "-------" "--------"
for EVENT in "${EVENTS[@]}"; do
    Q1=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q1_more_events.txt"       2>/dev/null || echo 0)
    Q2=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q2_zero_to_events_genes.txt" 2>/dev/null || echo 0)
    Q3=$(wc -l < "${OUTPUT_DIR}/${EVENT}_Q3_novel_transcripts.ioe"   2>/dev/null || echo 0)
    printf "  %-6s  %10s  %10s  %10s\n" "$EVENT" "$Q1" "$Q2" "$Q3"
done

echo ""
echo "  Output files saved to: ${OUTPUT_DIR}/"
echo "  ├── {EVENT}_Q1_more_events.txt         gene | strict | var1 | gained"
echo "  ├── {EVENT}_Q2_zero_to_events_genes.txt  gene_id list"
echo "  ├── {EVENT}_Q2_zero_to_events_counts.txt gene | var1_events"
echo "  └── {EVENT}_Q3_novel_transcripts.ioe    full IOE rows with UUID IDs"
echo ""
