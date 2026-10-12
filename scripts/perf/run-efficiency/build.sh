#!/bin/bash
# run-speed: build src/ into build/ (Release, FINETUNE=ON, stock OpenMS 3.5.0 env, the b1-ms1mem recipe),
# copy the binary to bin/DIALibGen-$TAG, then (CTEST=1) full CTest and the parameter-reference check.
# Never under the timing lock; nice 5.
set -e
export LC_ALL=C
A=/scratch/kohlbach/run-eff
TAG=${1:?tag}
ENV=/ceph/ibmi/abi/oliver/envs/dialibrefine
CXX=/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/deps/bin/x86_64-conda-linux-gnu-g++
MODELS=/ceph/ibmi/abi/oliver/opt/dialibgen-0.10.1/share/DIALibGen/models
PY=/ceph/ibmi/abi/oliver/envs/dialibgen-release-tests-0.11.0/bin/python
S=$A/src; B=$A/build
echo "uptime_before $(uptime)"
if [ ! -f $B/build.ninja ]; then
echo "configure $(date -u +%FT%TZ)"
$ENV/bin/cmake -G Ninja -S $S -B $B -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_MAKE_PROGRAM=$ENV/bin/ninja -DCMAKE_PREFIX_PATH=$ENV -DCMAKE_CXX_COMPILER=$CXX \
  -DTorch_DIR=/ceph/ibmi/abi/oliver/opt/libtorch-2.10-cpu-dialibgen/share/cmake/Torch \
  -DDIALIBGEN_BUILD_FINETUNE=ON -DODIA_MODEL_DIR=$MODELS -DDIALIBGEN_PEPTDEEP_ONNX_DIR=$MODELS \
  -DODIA_PYTHON=$PY -DDIALIBGEN_REFINEMENT_PYTHON=$PY -DDIALIBGEN_EXAMPLE_FASTA:FILEPATH=- > $A/logs/configure-$TAG.log 2>&1
fi
echo "build $(date -u +%FT%TZ)"
nice -n 19 $ENV/bin/ninja -C $B -j 8 > $A/logs/build-$TAG.log 2>&1 || { tail -60 $A/logs/build-$TAG.log; echo BUILD_FAILED; exit 1; }
echo "BUILD_OK $(date -u +%FT%TZ)"
cp $B/DIALibGen $A/bin/DIALibGen-$TAG && md5sum $A/bin/DIALibGen-$TAG
if [ "${CTEST:-0}" = 1 ]; then
  cd $B
  set +e
  # generate_under_thread_limit allows 12 threads above the user's current count: alone, not beside 7 other tests.
  nice -n 5 $ENV/bin/ctest -j 8 -E generate_under_thread_limit --output-on-failure > $A/logs/ctest-$TAG.log 2>&1
  echo "ctest rc=$?"
  nice -n 5 $ENV/bin/ctest -R generate_under_thread_limit --output-on-failure >> $A/logs/ctest-$TAG.log 2>&1
  echo "ctest thread-limit rc=$?"
  tail -12 $A/logs/ctest-$TAG.log
  $PY $S/scripts/update-parameters.py --check $B/DIALibGen > $A/logs/params-$TAG.log 2>&1; echo "params rc=$?"
  cat $A/logs/params-$TAG.log | tail -5
fi
echo "uptime_after $(uptime)"
echo "DONE_BUILD $TAG $(date -u +%FT%TZ)"
