#!/bin/bash
# gdbsample.sh LABEL INTERVAL_S THREADS "ENV" [FASTA]: UNTIMED diagnostic. gdb launches the knob build (traced, sub10 by
# default) and stops it every INTERVAL_S to record all threads' backtraces (diag/LABEL/s_<epoch>.txt; cost per stop in
# samples.tsv). The stops perturb the run, so its timings are not used. Run under: flock /scratch/kohlbach/bid10.lock
export LC_ALL=C
export DIALIBGEN_MODEL_DIR=/ceph/ibmi/abi/oliver/opt/dialibgen-0.10.1/share/DIALibGen/models
A=/scratch/kohlbach/bid10-threads; I=/scratch/kohlbach/bid10-instrument/in
L=$1; INT=$2; T=$3; E=$4; F=${5:-sub10}
D=$A/diag/$L; rm -rf $D; mkdir -p $D
{ echo "uptime_before $(uptime)"; echo "env $E threads $T fasta $F interval $INT"; echo "start $(date -u +%FT%TZ)"; } > $D/meta.txt
( while [ ! -s $D/pid ]; do sleep 0.2; done; P=$(cat $D/pid)
  while kill -0 $P 2>/dev/null; do sleep $INT; kill -INT $P 2>/dev/null; done ) &
K=$!
env $E ODIA_TRACE=$D/trace GDBS_DIR=$D GDBS_INT=$INT gdb -q -batch -x $A/tools/gdbsample.py --args \
  $A/bin/DIALibGen-knobs-ba7d248 -in $I/$F.fasta -config $I/G0.json -out $A/out/$L.parquet -threads $T > $D/out.log 2>&1
echo "rc $?" >> $D/meta.txt; kill $K 2>/dev/null; echo "end $(date -u +%FT%TZ)" >> $D/meta.txt
echo "uptime_after $(uptime)" >> $D/meta.txt
md5sum $A/out/$L.parquet >> $D/meta.txt; rm -f $A/out/$L.parquet
echo GDB_DONE $L $(date -u +%FT%TZ)
