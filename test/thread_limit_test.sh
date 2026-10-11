#!/usr/bin/env bash
# Under a thread limit generation must still finish, with the same library:
# a thread that cannot be started leaves the work to the threads that exist
# (parallel::startThreads). Until 0.11.2 the process aborted with
# "terminate called without an active exception".
#   thread_limit_test.sh DIALibGen FASTA WORKDIR      (Linux only; skips with 77 elsewhere)
set -u
tool=$1 fasta=$2 work=$3
[ "$(uname -s)" = Linux ] || { echo "not Linux; skipped"; exit 77; }
rm -rf "$work"; mkdir -p "$work"
# OpenBLAS starts its own threads before main and fails loudly under the limit.
export OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=
"$tool" -in "$fasta" -out "$work/free.tsv" -threads 64 > "$work/free.log" 2>&1 \
  || { cat "$work/free.log"; exit 1; }
# RLIMIT_NPROC counts every thread of the user on this machine: allow a dozen more.
running=$(ps -L -u "$(id -u)" --no-headers | wc -l)
limit=$((running + 12))
( ulimit -u "$limit" && exec "$tool" -in "$fasta" -out "$work/limited.tsv" -threads 64 ) > "$work/limited.log" 2>&1
rc=$?
if [ $rc -ne 0 ]; then
  echo "exit $rc under ulimit -u $limit (user threads: $running)"; tail -5 "$work/limited.log"; exit 1
fi
cmp "$work/free.tsv" "$work/limited.tsv" || { echo "the library differs under the thread limit"; exit 1; }
echo "ulimit -u $limit: same library"
