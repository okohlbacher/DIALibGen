#!/usr/bin/env python3
"""eqcheck.py REF.parquet OTHER.parquet [...]: is each OTHER's table equal to REF's (pyarrow Table.equals, data and
schema without metadata), and which schema-metadata entries differ (JSON values compared with 'tool_version' removed)."""
import sys, json
import pyarrow.parquet as pq

def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k != 'tool_version'}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o

def meta(t):
    out = {}
    for k, v in (t.schema.metadata or {}).items():
        k = k.decode(); v = v.decode()
        try:
            out[k] = strip(json.loads(v))
        except Exception:
            out[k] = v
    return out

ref = pq.read_table(sys.argv[1])
mr = meta(ref)
for p in sys.argv[2:]:
    t = pq.read_table(p)
    eq = t.equals(ref)
    m = meta(t)
    diff = sorted(k for k in set(m) | set(mr) if m.get(k) != mr.get(k))
    print('%s\trows=%d\tTable.equals=%s\tmetadata_diff(after tool_version removed)=%s' % (p, t.num_rows, eq, ','.join(diff) or 'none'))
