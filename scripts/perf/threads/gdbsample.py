# gdb -batch -x gdbsample.py --args PROG ...   (env GDBS_DIR, GDBS_INT)
# kernel.yama.ptrace_scope=1 on kim: gdb must be the program's ancestor, so gdb launches it and a timer thread
# is interrupted with SIGINT (not passed to the program) every GDBS_INT s by gdbsample.sh; at each stop every thread's backtrace is written.
import gdb, os, threading, time, signal
D = os.environ['GDBS_DIR']; INT = float(os.environ.get('GDBS_INT', '2'))
for c in ("set pagination off", "set confirm off", "set print thread-events off",
          "handle SIGINT stop print nopass", "handle SIGPIPE nostop noprint pass", "set width 0"):
    gdb.execute(c)
log = open(os.path.join(D, 'samples.tsv'), 'w'); log.write('epoch_s\tcost_s\tthreads\n')
def alive():
    inf = gdb.selected_inferior()
    return inf is not None and inf.pid != 0
gdb.execute("starti", to_string=True)
pid = gdb.selected_inferior().pid
# gdb holds the GIL inside gdb.execute, so the periodic SIGINT comes from an outside shell loop (gdbsample.sh)
open(os.path.join(D, 'pid'), 'w').write(str(pid))
try:
    gdb.execute("continue", to_string=True)
except gdb.error:
    pass
while alive():
    t0 = time.time()
    try:
        out = gdb.execute("thread apply all bt 14", to_string=True)
    except gdb.error as e:
        out = 'ERR ' + str(e)
    with open(os.path.join(D, 's_%.3f.txt' % t0), 'w') as f: f.write(out)
    log.write('%.3f\t%.3f\t%d\n' % (t0, time.time() - t0, out.count('\nThread ') + out.startswith('Thread ')))
    log.flush()
    try:
        gdb.execute("continue", to_string=True)
    except gdb.error:
        break
log.close()
