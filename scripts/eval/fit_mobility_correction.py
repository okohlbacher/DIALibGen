"""Fit the straight-line correction of DIALibGen's derived 1/K0 (Mason-Schamp from the predicted CCS) to the observed 1/K0
of diaPASEF runs: the mobility_correction lines of data/instruments/timstof.yaml.

LIB.parquet is a DIALibGen library made WITHOUT the correction (-generation:mobility_correction off; 0.12.1's timsTOF
library); each report.parquet a DIA-NN 2.x search of one training run (precursors at Q.Value <= 0.01, targets, charges
2-4). RT_IM_VALUES.parquet is an evaluation reference (columns run, kind, library, mode, pid, heldout, obs, pred: the
timsTOF benchmark's rt_im.py output) on which the lines are scored, rows kind=im, library=G0, mode=stored.

Models: global  1/K0_obs = a + b * 1/K0_pred
        charge  the same per precursor charge (2, 3, 4)
Fitted by least squares after trimming residuals beyond 4 robust SD (two passes); every training run weighs the same.
usage: fit_mobility_correction.py LIB.parquet RT_IM_VALUES.parquet OUT_PREFIX name=report.parquet ...
Writes OUT_PREFIX.json (the lines) and OUT_PREFIX.eval.tsv (median |d1/K0|, bias, by charge and 1/K0 bin)."""
import sys, json
import numpy as np, pandas as pd, pyarrow.parquet as pq

lib_path, ref_path, out = sys.argv[1:4]
train = [a.split("=", 1) for a in sys.argv[4:]]


def norm(s):
    return s.str.replace("(Carbamidomethyl)", "(UniMod:4)", regex=False).str.replace("(Oxidation)", "(UniMod:35)", regex=False)


lib = pq.read_table(lib_path, columns=["Precursor.Id", "IM", "Precursor.Charge"]).to_pandas().drop_duplicates("Precursor.Id")
lib.index = norm(lib["Precursor.Id"])
rows = []
for name, path in train:
    r = pq.read_table(path, columns=["Precursor.Id", "IM", "Q.Value", "Decoy", "Precursor.Charge"]).to_pandas()
    r = r[(r["Q.Value"] <= 0.01) & (r.Decoy == 0) & (r["Precursor.Charge"].between(2, 4))].drop_duplicates("Precursor.Id")
    r = r[r["Precursor.Id"].isin(lib.index)]
    rows.append(pd.DataFrame({"run": name, "z": r["Precursor.Charge"].values, "obs": r.IM.values,
                              "pred": lib.loc[r["Precursor.Id"], "IM"].values}))
T = pd.concat(rows, ignore_index=True)
T = T[np.isfinite(T.obs) & np.isfinite(T.pred) & (T.obs > 0)]


def fit_line(d):
    keep = np.ones(len(d), bool)
    w = d.run.map(1.0 / d.run.value_counts()).values
    for _ in range(3):
        W = np.sqrt(w[keep])
        A = np.column_stack([np.ones(keep.sum()), d.pred.values[keep]]) * W[:, None]
        a, b = np.linalg.lstsq(A, d.obs.values[keep] * W, rcond=None)[0]
        res = d.obs.values - a - b * d.pred.values
        sd = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        keep = np.abs(res) <= 4 * sd
    return float(a), float(b), int(keep.sum())


models = {"none": {}, "global": {}, "charge": {}}
a, b, n = fit_line(T); models["global"] = {z: (a, b) for z in (2, 3, 4)}; fitinfo = {"global": dict(a=a, b=b, n=n)}
for z in (2, 3, 4):
    a, b, n = fit_line(T[T.z == z]); models["charge"][z] = (a, b); fitinfo[f"charge{z}"] = dict(a=a, b=b, n=n)
models["none"] = {z: (0.0, 1.0) for z in (2, 3, 4)}
print("training:", T.groupby(["run", "z"]).size().to_dict())
print("fits:", json.dumps(fitinfo, indent=1))
# training diagnostics: median residual by run and charge before/after
for m in ("none", "global", "charge"):
    c = np.array([models[m][z][0] for z in T.z]) + np.array([models[m][z][1] for z in T.z]) * T.pred.values
    d = c - T.obs.values
    print(f"train {m:7s} median d {np.median(d):+.4f} median|d| {np.median(np.abs(d)):.4f}",
          {f"{run}/{z}": round(float(np.median(d[(T.run == run).values & (T.z == z).values])), 4) for run in T.run.unique() for z in (2, 3, 4)})

# evaluation on the benchmark reference
R = pq.read_table(ref_path).to_pandas()
R = R[(R.kind == "im") & (R.library == "G0") & (R["mode"] == "stored")].copy()
R["z"] = R.pid.str.extract(r"(\d+)$")[0].astype(int)
res = []
for m in ("none", "global", "charge"):
    corr = np.array([models[m][z][0] for z in R.z]) + np.array([models[m][z][1] for z in R.z]) * R.pred.values
    R[m] = corr - R.obs.values
    for sel, mask in (("run01 held-out", (R.run == 1) & R.heldout), ("runs 02+03", R.run > 1)):
        d = R.loc[mask, m].values
        row = dict(model=m, set=sel, n=len(d), median_abs=np.median(np.abs(d)), bias=np.median(d),
                   p95=np.percentile(np.abs(d), 95), within_001=(np.abs(d) <= 0.01).mean())
        for z in (2, 3, 4):
            row[f"bias_z{z}"] = np.median(R.loc[mask & (R.z == z), m])
        bins = pd.cut(R.loc[mask, "obs"], [0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.6])
        for k, v in R.loc[mask, m].groupby(bins, observed=True).median().items():
            row[f"bias_{k.left:.1f}-{k.right:.1f}"] = v
        res.append(row)
E = pd.DataFrame(res)
pd.set_option("display.width", 250)
print(E.round(4).to_string())
E.to_csv(out + ".eval.tsv", sep="\t", index=False, float_format="%.5f")
json.dump({"fits": fitinfo, "training_rows": int(len(T))}, open(out + ".json", "w"), indent=1)
