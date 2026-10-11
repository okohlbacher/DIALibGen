#!/bin/bash
# smoke.sh BIN TAG -- untimed correctness smoke on tiny.fasta, 4 threads, traced; run under the lock
export LC_ALL=C
export DIALIBGEN_MODEL_DIR=/ceph/ibmi/abi/oliver/opt/dialibgen-0.10.1/share/DIALibGen/models
A=/scratch/kohlbach/bid10-solutions; IN=/scratch/kohlbach/bid10-instrument/in; BIN=$1; T=${2:-v1}
mkdir -p $A/smoke; cd $A/smoke; rm -rf $T; mkdir $T
for v in base:"" sol:"ODIA_SOLUTIONS=1" solprune:"ODIA_SOLUTIONS=1 ODIA_MS2_PRUNED=$A/models/peptdeep_ms2_pruned.onnx"; do
  n=${v%%:*}; e=${v#*:}
  for ext in parquet tsv; do
    env ODIA_TRACE=$A/smoke/$T/tr_${n}_$ext $e $BIN -in $IN/tiny.fasta -config $IN/G0.json -out $A/smoke/$T/${n}.$ext -threads 4 > $A/smoke/$T/${n}_$ext.log 2>&1
    echo "$n $ext rc=$? $(md5sum < $A/smoke/$T/${n}.$ext | cut -c1-32)"
  done
done
echo SMOKE_DONE
