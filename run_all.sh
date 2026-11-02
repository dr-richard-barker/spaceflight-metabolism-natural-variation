#!/usr/bin/env bash
# Regenerate every table, figure and manuscript number from the vendored inputs.
#
#   ./run_all.sh            analysis + numbers.tex
#   ./run_all.sh --vendor   re-vendor inputs from the sibling repositories first
#   ./run_all.sh --pdf      also build manuscript/latex/main.pdf and main.docx
#
# Not re-run here, on purpose:
#   scripts/11_prereg.py,           the pre-registrations are frozen; re-running them would rewrite
#   scripts/21_prereg_osd38.py      prereg/ from today's KEGG release
#   scripts/14_submit_paintomics.py PaintOmics jobs are submitted once; their records are vendored
#                                   in results/paintomics_raw/ and read from there
#   scripts/27_... versions          the OSDR file listing is a dated snapshot
#                                   (data/reference/osd38_osdr_versions.tsv)
# Revision 1 analyses (scripts/23-28, prereg/revision1) also need R with limma and edgeR (Rscript on PATH).
# Network: the Wu supplement, KEGG lists, MetaboLights, OSDR (OSD-38, OSD-136, OSD-678) and GEO (GSE14502,
# GPL198) files are fetched into the git-ignored data/cache/ on first run and checked against their
# recorded SHA-256 (data/reference/*_sources.tsv).
set -euo pipefail
cd "$(dirname "$0")/scripts"

for arg in "$@"; do
  if [ "$arg" = "--vendor" ]; then python3 00_vendor_inputs.py; python3 20_osd38_inputs.py contrasts; fi
done

python3 10_light_dark_reference.py fetch
python3 10_light_dark_reference.py build
python3 12_model_vs_reference.py
python3 20_osd38_inputs.py fetch
python3 22_osd38_tests.py
python3 13_build_layer.py
python3 15_paintomics_tables.py
python3 16_gwas_light_dark_tests.py
python3 17_gcms_power.py fetch
python3 17_gcms_power.py build
python3 23_reference_strata.py
python3 24_osd678_interaction.py
python3 25_geneset_diagnostics.py
python3 26_hypoxia_positive_control.py
python3 27_osd38_proteome_checks.py build
python3 28_paintomics_controls.py
python3 18_manuscript_numbers.py
python3 make_manifest.py

for arg in "$@"; do
  if [ "$arg" = "--pdf" ]; then
    (cd ../manuscript/latex && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex && ./build_docx.sh)
  fi
done
python3 19_build_site.py
echo "run_all: done"
