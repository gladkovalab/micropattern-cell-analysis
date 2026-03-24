# 20251229_paths_for_analysis.txt was generated using the following command
# find /groups/vale/valelab/_for_Mark/patterned_data/ -maxdepth 2 -mindepth 2 -type d | grep -v example > 20251229_paths_for_analysis.txt
for _dir in `cat 20251229_paths_for_analysis.txt`
do
    echo $_dir
    bsub -n 8 -P vale pixi run python template_matching_bulk.py $_dir
done
