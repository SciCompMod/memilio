#!/bin/bash
# Runs all arms of the location split study one after another (finished arms are skipped), then the
# paired analysis and the figures. All arms use the calibrated defaults of the simulation unless noted.
D=${STUDY_DIR:-$PWD/study_output}; B=${SIM_BINARY:-$PWD/build/bin/location_split_simulation}
PY=${PYTHON:-python3}
HERE=$(cd "$(dirname "$0")" && pwd); export STUDY_DIR=$D
OUT=${STUDY_OUT:-$D/results}; export STUDY_OUT=$OUT; mkdir -p $D $OUT

run() { # name runs options...
  local name=$1 runs=$2; shift 2
  grep -q "^$name " $D/progress.txt 2>/dev/null && return
  rm -rf $D/$name
  $B --runs $runs --days 90 --n_persons 10000 --initial_infections 10 --seed 42 ${THREADS:+--threads $THREADS} "$@" \
     --output_dir $D/$name > $D/$name.log 2>&1 || { echo "$name failed, see $D/$name.log"; exit 1; }
  echo "$name done $(date +%T)" >> $D/progress.txt; cp $D/progress.txt $OUT/progress.txt
  $PY $HERE/analyze.py 2>> $D/analyze.err
}

# School split at the default school weight: fixed rooms, and rooms redrawn every day. All school arms
# also write the per room and hour counts for the comparison with the one-step theory (theory_check.py).
for rs in 0 50 25 10 5; do
  run sw1_$([ $rs -eq 0 ] && echo no_split || echo rooms_$rs) 50 --school_room_size $rs --write_school_rooms
done
for rs in 50 25 10 5; do
  run sw1_redrawn_rooms_$rs 50 --school_room_size $rs --redraw_school_rooms --write_school_rooms
done

# Workplace and social event splits, only if WORK_EVENT_SPLITS=1.
if [ "${WORK_EVENT_SPLITS:-0}" = 1 ]; then
  for rs in 5 3; do
    run work_rooms_$rs 30 --work_room_size $rs
    run event_rooms_$rs 30 --event_room_size $rs
  done
fi

$PY $HERE/theory_check.py >> $D/theory_check.log 2>&1
$PY $HERE/make_figures.py >> $D/figures.log 2>&1
echo "pipeline finished $(date +%T)" >> $D/progress.txt; cp $D/progress.txt $OUT/progress.txt
