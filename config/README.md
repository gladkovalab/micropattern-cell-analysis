# config/

Configuration and input data files for the analysis pipeline.

## Files

### `Comparisons_table_v3.xlsx`
Defines experimental conditions and well assignments across plates, used by
`generate_comparison_projections.py` and `prototype_comparisons.py` via `comparison_loader.py`.

Original location: `/groups/vale/valelab/_for_Mark/analysis/Comparisons_table_v3.xlsx`

### `20251229_paths_for_analysis.txt`
List of well directories submitted for batch analysis via `bsub_analysis.sh`.
Generated with:
```bash
find /groups/vale/valelab/_for_Mark/patterned_data/ -maxdepth 2 -mindepth 2 -type d | grep -v example > config/20251229_paths_for_analysis.txt
```
