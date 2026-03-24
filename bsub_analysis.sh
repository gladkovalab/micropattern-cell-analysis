#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C03_250718_TRAK1/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C05_250718_TRAK1_mDRH/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250710_patterned_plate_9_good/C06_250718_TRAK1_mDRH_dSp/
  
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E05_250808_TRAK1_wt/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E06_250808_TRAK1_mDRH/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/E07_250811_TRAK2_mDRH_Sp/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250731_patterned_plate_11_good/F06_250811_TRAK1_mDRH_dSp/

#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/B06_TRAK1_wt_combined/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/E06_250606_TRAK1-mDRH/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250521_patterned_plate_1/D07_250606_TRAK1-mDRH+dSpindly/

#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B03_TRAK1_250616/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B06_250617_TRAK1_mDRH_dSp/
#bsub -n 8 -P vale pixi run python template_matching_bulk.py /groups/vale/valelab/_for_Mark/patterned_data/250612_patterned_plate_3/B05_250617_TRAK1_mDRH/

# find /groups/vale/valelab/_for_Mark/patterned_data/ -maxdepth 2 -mindepth 2 -type d | grep -v example > 20251229_paths_for_analysis.txt
for _dir in `cat 20251229_paths_for_analysis.txt`
do
    echo $_dir
    bsub -n 8 -P vale pixi run python template_matching_bulk.py $_dir
done
