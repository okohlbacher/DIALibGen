#!/bin/bash
# queue2.sh FILE [WAITPID]: like queue.sh, after process WAITPID (the previous queue) has exited.
# A line "CMD <shell command>" runs that command (not under the timing lock: builds, CTest);
# any other line is "VAR=.. [VAR=..] CELL TAG BIN THREADS DATA [args]" for cell.sh.
A=/scratch/kohlbach/run-eff
Q=$(readlink -f $1); W=${2:-}
if [ -n "$W" ]; then while [ -d /proc/$W ]; do sleep 15; done; fi
while true; do
  L=$(head -1 $Q 2>/dev/null)
  [ -z "$L" ] && break
  sed -i 1d $Q
  echo "$(date -u +%FT%TZ) RUN $L" >> $A/logs/queue.log
  case "$L" in
    CMD\ *) ( eval "${L#CMD }" ) >> $A/logs/queue.log 2>&1 ;;
    *) eval "env ${L%% CELL *} bash $A/tools/cell.sh ${L#* CELL }" >> $A/logs/queue.log 2>&1 ;;
  esac
done
echo "$(date -u +%FT%TZ) QUEUE_EMPTY $Q" >> $A/logs/queue.log
