# boundary_analysis

Author: Gricey

## What this is

This is a side investigation I ran into SUPPA's **boundary parameter** for
`generateEvents` — the setting that controls how strictly a splice site has
to match before SUPPA calls it the same alternative splicing event. `strict`
mode requires an exact coordinate match; `variable Nnt` mode tolerates a
shift of up to N nucleotides at the splice site.

The question I was answering: how much does relaxing this boundary change
what SUPPA detects, and does it change the biology we'd report? Kiran and
Khushi wanted to know before we committed to one setting for the actual
splicing analyses, since a boundary that's too strict could be silently
dropping real, biologically relevant events (especially retained introns),
while one that's too loose could inflate the event catalogue with
coordinate noise.

**Dataset:** samples KM1–KM9, mouse liver, 3 conditions x 3 replicates —
`CT8` (KM1–3), `CT20` (KM4–6), `PerKO` (KM7–9). Quantified with Salmon
against `transcriptome_productivity.gtf` / the same reference transcriptome
Khushi built for the PERIOD manuscript work. That shared reference is the
**only** thing this project has in common with the PERIOD manuscript
(`orthogonal_validation/`) — the KM1-9 dataset itself is separate and isn't
a source for manuscript figures.

## Workflow

1. **`scripts/run_suppa.sh`** — runs the full SUPPA pipeline (TPM prep →
   `generateEvents` → `psiPerEvent` → `diffSplice` → `clusterEvents`) for one
   boundary mode at a time. I reran it with `BOUNDARY=S` (strict) and
   `BOUNDARY=V` with `THRESHOLD=1,3,5,10,20` (variable Nnt), producing the
   `data/output_strict/` and `data/output_variable_*nt/` folders.

2. **`notebooks/01_filter_events.ipynb`** — for strict vs variable 1nt,
   works out which events are shared between modes vs newly detected by
   relaxing the boundary, at both the event level and the gene level.
   Includes a superset check (every strict event must also exist in
   variable) as a correctness proof of the coordinate-translation logic.

3. **`notebooks/02_barplots.ipynb`** — visualizes event counts across all
   six boundary modes (strict, 1/3/5/10/20nt) for all seven SUPPA event
   types, and shows that counts plateau after ~1nt — which is why variable
   1nt was chosen as the working parameter rather than something looser.

4. **`notebooks/03_psi_comparison.ipynb`** — for events detected in *both*
   modes, computes delta PSI (`|PSI_strict - PSI_variable|`) per event and
   plots its density. Also compares significant-event counts (CT20 vs
   PerKO) between modes, and picks the single highest-delta-PSI event as an
   IGV candidate.

5. **`scripts/compare_event_populations.sh`**, **`compare_runs.sh`**,
   **`get_igv_candidates.sh`**, **`igv_filter.sh`** — command-line helpers
   used alongside the notebooks: tallying gene/event populations across
   modes, comparing significant-event counts across all six boundary runs
   at once, and two different ways of shortlisting IGV candidates (by
   dPSI/p-value significance, or by a clean small jump in event count per
   gene).

6. **`figures/igv/`** — manual IGV screenshots of the candidate events
   picked above, used to visually confirm that the extra events detected in
   variable mode are real alternative splicing and not coordinate noise.

## Folder layout

- `data/output_strict/`, `data/output_variable_{1,3,5,10,20}nt/` — one full
  SUPPA run per boundary mode (`events/`, `tpm/`, `psi/`, `diff/`, `cluster/`)
- `data/ioe_populations/`, `data/ioe_diff/` — Q1–Q3 gene/event population
  comparisons and the per-sample delta PSI table (`dpsi_table_KM4.csv`)
- `data/salmon/` — per-sample Salmon quant output (KM1–KM9)
- `figures/plots/` — notebook-generated bar plots, density plots, Venn
  diagrams, summary tables
- `figures/igv/Tests/` and `figures/igv/Comparisons/` — IGV screenshots,
  early exploratory vs the final candidate comparisons
- `notebooks/` — the analysis notebooks described above
- `scripts/` — the pipeline and comparison scripts described above

## Next step

The density plot of delta PSI (`figures/plots/density_dpsi.png`, notebook
03 Step 3) has a peak near 0 (strict and variable agree) plus a smaller
population out in the tails, where the two boundary modes disagree a lot on
PSI for the same event. I haven't looked at those tail genes specifically
yet — next thing to do is pull out the events/genes sitting in the tails of
that distribution and quantify them properly (how many, which genes, how
large is "large"), rather than just eyeballing the single highest-delta
example I picked for IGV in Step 6.
