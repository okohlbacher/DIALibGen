#!/bin/bash
# cell.sh TAG BIN THREADS DATA [extra DIALibGen args]: one timed, traced -run under the dax timing lock.
#   DATA tims       the reduced timsTOF cell of b1-ms1mem/runs/cell_auto (30,000 pairs, library intensities, -mode refine)
#        tims_pred  the same with search:intensities predicted (both members predicted by the MS2 model)
#        astral     reduced Astral cell: full Astral run, 40,000 pairs, library intensities, -mode refine
#        full       the timsTOF -mode tune command of b1-ms1mem/runs/tims_auto
#   env  CAP=N      ODIA_OPENSWATH_THREADS=N (the OpenSWATH team cap)
#        OMPENV="K=V ..."  exported before launch (e.g. OMP_WAIT_POLICY=passive GOMP_SPINCOUNT=10000)
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
  full) OUT=$D/tuned.tsv
        ARGS=(-mode tune -tune_heads both -search:entrapment_tag ENTRAP_ -in $GLIB -run $TIMS -out $OUT -threads $T
              -search:ms1_in_memory auto) ;;
  tune) OUT=$D/tuned.tsv
        ARGS=(-mode tune -tune_heads both -in $GLIB -ids /scratch/kohlbach/run-speed/runs/final_full_t16/tuned.tsv.ids.parquet
              -out $OUT -threads $T) ;;
  *) echo "unknown DATA $DATA"; exit 2 ;;
esac
ARGS+=("$@")
export ODIA_TRACE=$D/trace
[ -n "${CAP:-}" ] && export ODIA_OPENSWATH_THREADS=$CAP
for kv in ${OMPENV:-}; do export "$kv"; done
{ echo "# $(date -u +%FT%TZ) data=$DATA threads=$T bin_md5=$(md5sum $BIN | cut -c1-32)"
  env | grep -E '^(ODIA_|OMP_|GOMP_|OPENBLAS_|DIALIBGEN_MODEL_DIR)' | sort | sed 's/^/# env /'
  echo "PIN=${PIN:-} $BIN ${ARGS[*]}"; } > cmd.txt
exec 9>/scratch/kohlbach/dax-timing.lock
echo "waiting for lock $(date -u +%FT%TZ)" > cond.txt
flock 9
{ echo "start $(date -u +%FT%TZ) $(date +%s)"; uptime; free -g | head -2; } >> cond.txt
/usr/bin/time -v -o time.txt ${PIN:-} nice -n 5 $BIN "${ARGS[@]}" > log.txt 2>&1
echo "rc=$?" >> log.txt
{ echo "end $(date -u +%FT%TZ) $(date +%s)"; uptime; } >> cond.txt
flock -u 9
exec 9>&-
grep -E "Elapsed|Maximum resident|Exit status" time.txt >> log.txt
rm -rf $D/.dialibgen-run-cache-* 2>/dev/null
md5sum $OUT.ids.parquet $OUT > md5.txt 2>/dev/null
[ "$DATA" = full ] && rm -f $OUT
echo CELL_DONE $TAG $(date -u +%FT%TZ)
