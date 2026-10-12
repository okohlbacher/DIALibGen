#!/bin/bash
# gdbcell.sh TAG BIN THREADS DATA [args]: a traced DIAGNOSTIC -run (not used for timing) under the dax timing lock,
# with gdb stack samples of all threads (LD_PRELOAD ptracer.so allows the attach under ptrace_scope 1).
#   env START=regex  start sampling once log.txt matches (default: at once)
#       STOP=regex   stop sampling once log.txt matches; KILL=1 then kills the run (default: sample to the end)
#       GAP=s        seconds between samples (default 10); NS=n at most n samples (default 200)
#       CAP=N, OMPENV as cell.sh
export LC_ALL=C OPENBLAS_NUM_THREADS=1
A=/scratch/kohlbach/run-eff
TAG=$1; BIN=$2; T=$3; DATA=$4; shift 4
TIMS=/scratch/kohlbach/mzpc-finish/clean/new-full.mzML
GLIB=/scratch/kohlbach/bid3-accept/tims/G0.pg.tsv
ASTRAL=/ceph/ibmi/abi/dont-backup/kohlbach/odia/data/astral.mzML
ALIB=/scratch/kohlbach/bid-accept/lib/library.parquet
export DIALIBGEN_MODEL_DIR=/ceph/ibmi/abi/oliver/opt/dialibgen-0.11.0/share/DIALibGen/models
D=$A/runs/$TAG
rm -rf $D && mkdir -p $D && cd $D
case $DATA in
  tims) OUT=$D/refined.parquet
        ARGS=(-mode refine -in $GLIB -run $TIMS -out $OUT -threads $T -search:entrapment_tag ENTRAP_ -search:max_pairs 30000
              -search:intensities library -search:cache_dir $D -search:ms1_in_memory auto) ;;
  tims_pred) OUT=$D/refined.parquet
        ARGS=(-mode refine -in $GLIB -run $TIMS -out $OUT -threads $T -search:entrapment_tag ENTRAP_ -search:max_pairs 30000
              -search:cache_dir $D -search:ms1_in_memory auto) ;;
  astral) OUT=$D/refined.parquet
        ARGS=(-mode refine -in $ALIB -run $ASTRAL -out $OUT -threads $T -search:max_pairs 40000 -search:intensities library
              -search:cache_dir $D) ;;
  *) echo "unknown DATA $DATA"; exit 2 ;;
esac
ARGS+=("$@")
export ODIA_TRACE=$D/trace
[ -n "${CAP:-}" ] && export ODIA_OPENSWATH_THREADS=$CAP
for kv in ${OMPENV:-}; do export "$kv"; done
{ echo "# DIAGNOSTIC (gdb-sampled, not timed) $(date -u +%FT%TZ) data=$DATA threads=$T bin_md5=$(md5sum $BIN | cut -c1-32)"
  echo "# START=${START:-} STOP=${STOP:-} KILL=${KILL:-0} GAP=${GAP:-10} NS=${NS:-200}"
  env | grep -E '^(ODIA_|OMP_|GOMP_|OPENBLAS_|DIALIBGEN_MODEL_DIR)' | sort | sed 's/^/# env /'
  echo "$BIN ${ARGS[*]}"; } > cmd.txt
exec 9>/scratch/kohlbach/dax-timing.lock
echo "waiting for lock $(date -u +%FT%TZ)" > cond.txt
flock 9
{ echo "start $(date -u +%FT%TZ) $(date +%s)"; uptime; free -g | head -2; } >> cond.txt
LD_PRELOAD=$A/tools/ptracer.so nice -n 5 $BIN "${ARGS[@]}" > log.txt 2>&1 &
P=$!
if [ -n "${START:-}" ]; then until grep -qE "$START" log.txt 2>/dev/null; do sleep 0.5; kill -0 $P 2>/dev/null || break; done; fi
k=0
while kill -0 $P 2>/dev/null && [ $k -lt ${NS:-200} ]; do
  if [ -n "${STOP:-}" ] && grep -qE "$STOP" log.txt 2>/dev/null; then break; fi
  k=$((k+1))
  echo "=== sample $k $(date +%s.%N)" >> stacks.txt
  timeout 60 gdb -q -batch -p $P -ex "set pagination off" -ex "thread apply all bt 24" >> stacks.txt 2>/dev/null
  echo "=== end $k $(date +%s.%N)" >> stacks.txt
  sleep ${GAP:-10}
done
[ "${KILL:-0}" = 1 ] && kill $P 2>/dev/null
wait $P; echo "rc=$?" >> log.txt
{ echo "end $(date -u +%FT%TZ) $(date +%s)"; uptime; } >> cond.txt
flock -u 9
exec 9>&-
rm -rf $D/.dialibgen-run-cache-* 2>/dev/null
md5sum $OUT.ids.parquet $OUT > md5.txt 2>/dev/null
echo GDB_DONE $TAG $k samples $(date -u +%FT%TZ)
