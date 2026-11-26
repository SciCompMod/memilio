#!/bin/bash
# Sequential pipeline: main school arms, work/event splits, school-share sensitivity, figures, presentation.
D=${STUDY_DIR:-$PWD/study_output}; B=${SIM_BINARY:-$PWD/build/bin/location_split_simulation}
HERE=$(cd "$(dirname "$0")" && pwd); export STUDY_DIR=$D
OUT=${STUDY_OUT:-$D/results}; export STUDY_OUT=$OUT; mkdir -p $D $OUT
W="home=20.3425,work=0.625,social_event=6.08233,shop=0.1584"
run() { # name runs options...
  local name=$1 runs=$2; shift 2
  grep -q "^$name " $D/progress.txt 2>/dev/null && return
  rm -rf $D/$name
  $B --runs $runs --days 90 --n_persons 10000 --initial_infections 10 --seed 42 ${THREADS:+--threads $THREADS} "$@" --output_dir $D/$name > $D/$name.log 2>&1
  echo "$name done $(date +%T)" >> $D/progress.txt; cp $D/progress.txt $OUT/progress.txt
  python3 $HERE/analyze.py 2>> $D/analyze.err
}
run sw1_rooms_25 50 --infection_k 1.596 --contact_weights $W,school=1.6 --school_room_size 25
run sw1_rooms_50 50 --infection_k 1.596 --contact_weights $W,school=1.6 --school_room_size 50
run work_rooms_5 30 --infection_k 1.596 --contact_weights $W,school=1.6 --work_room_size 5
run event_rooms_5 30 --infection_k 1.596 --contact_weights $W,school=1.6 --event_room_size 5
[ -f $D/scan/scan_result.json ] || python3 $HERE/scan_school.py >> $D/scan.out 2>&1
for t in 0.15 0.25; do
  tag=share$(python3 -c "print(int(round($t*100)))")
  read w k <<< $(python3 -c "import json; r=json.load(open('$D/scan/scan_result.json'))['targets']['$t']; print(f\"{r['school_weight']:.4f} {r['k']:.4f}\")")
  for rs in 0 25 5; do
    name=${tag}_$([ $rs -eq 0 ] && echo no_split || echo rooms_$rs)
    run $name 30 --infection_k $k --contact_weights $W,school=$w --school_room_size $rs
  done
done
run work_rooms_3 30 --infection_k 1.596 --contact_weights $W,school=1.6 --work_room_size 3
run event_rooms_3 30 --infection_k 1.596 --contact_weights $W,school=1.6 --event_room_size 3
python3 $HERE/make_figures.py >> $D/figures.log 2>&1
echo "pipeline finished $(date +%T)" >> $D/progress.txt; cp $D/progress.txt $OUT/progress.txt
