#!/bin/bash
# Stage the six missing TRAK2 wells from SMB to local disk using parallel cp.
# `cp` is generally faster than `rsync` for one-off bulk transfers since it
# doesn't compute checksums. We run 3 wells in parallel (SMB server can
# usually handle 3-4 concurrent streams without flapping).
#
# Local layout:
#   replication/local_staged/patterned_data/{plate}/{well_dir}/Cell*.nd2
set -u
SRC=/Volumes/valelab/_for_Mark/patterned_data
DST=/Users/gladkoc/Dev/micropattern_cell_analysis/replication/local_staged/patterned_data
LOG=/Users/gladkoc/Dev/micropattern_cell_analysis/replication/local_staged/stage.log
mkdir -p "$DST"
: > "$LOG"

wells=(
  "250606_patterned_plate_2/D05_250612_TRAK2-mDRH+mSpin"
  "250612_patterned_plate_3/B07_TRAK2_mDRH"
  "250612_patterned_plate_3/B08_TRAK2_mDRH_dSp"
  "250710_patterned_plate_9_good/C07_250718_TRAK2_mDRH"
  "250710_patterned_plate_9_good/C08_250721_TRAK2_mDRH_dSp"
  "250731_patterned_plate_11_good/D07_250811_TRAK2_mDRH"
  "250731_patterned_plate_11_good/E07_250811_TRAK2_mDRH_mSpin"
)

stage_one() {
  local rel="$1"
  local src_matches=( "$SRC"/"${rel%/*}"/"${rel##*/}"* )
  local src="${src_matches[0]}"
  if [ ! -d "$src" ]; then
    echo "[$(date +%H:%M:%S)] SKIP $rel (source dir not found)" >> "$LOG"
    return
  fi
  local dst_plate="$DST/${rel%/*}"
  mkdir -p "$dst_plate"
  local dst="$dst_plate/$(basename "$src")"
  mkdir -p "$dst"
  echo "[$(date +%H:%M:%S)] START $(basename "$src")" >> "$LOG"
  local ok=0 fail=0
  # Copy only top-level Cell*.nd2 / cell*.nd2; skip subdirs via shell globs.
  for f in "$src"/Cell*.nd2 "$src"/cell*.nd2; do
    [ -f "$f" ] || continue
    local base=$(basename "$f")
    local out="$dst/$base"
    if [ -f "$out" ] && [ "$(stat -f %z "$out" 2>/dev/null)" = "$(stat -f %z "$f" 2>/dev/null)" ]; then
      echo "[$(date +%H:%M:%S)]   SKIP $(basename "$src")/$base (size match)" >> "$LOG"
      ((ok++))
      continue
    fi
    echo "[$(date +%H:%M:%S)]   cp $(basename "$src")/$base …" >> "$LOG"
    if cp "$f" "$out.part" 2>>"$LOG" && mv "$out.part" "$out"; then
      ((ok++))
      echo "[$(date +%H:%M:%S)]   DONE $(basename "$src")/$base" >> "$LOG"
    else
      ((fail++))
      rm -f "$out.part"
      echo "[$(date +%H:%M:%S)]   FAIL $(basename "$src")/$base" >> "$LOG"
    fi
  done
  echo "[$(date +%H:%M:%S)] WELL DONE $(basename "$src")  ok=$ok fail=$fail" >> "$LOG"
}
export -f stage_one
export SRC DST LOG

# Parallel — 3 wells at a time
printf '%s\n' "${wells[@]}" | xargs -n1 -P3 -I{} bash -c 'stage_one "$@"' _ {}

echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
du -sh "$DST" | tee -a "$LOG"
