#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -f data/proteome_records.json || ! -f results/metadata.json ]]; then
  echo "Raw UniProt snapshot not found; downloading via UniProt REST API."
  python src/fetch_proteome.py
fi
python src/parse_annotations.py
python src/calculate_properties.py
python src/select_candidates.py
python src/check_reference_enzymes.py
python src/audit_pipeline.py
python src/make_figures.py
python -m pytest tests/ -q
printf '\nPipeline, report, plots and tests finished.\n'
