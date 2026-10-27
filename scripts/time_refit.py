"""Full-dev refit timing of the G4 finalists (LESSONS 14: time on a quiet machine).

    python scripts/time_refit.py --models m5_avg,m1_E_resB0 --threads 4

Each model is refitted on all 1022 DepMap dev lines with the G4 protocol (scripts/heldout_common.refit_predict:
scripts/run.py tuning, seeds, leakage sentinel) and predicts ONE DepMap noncancerous (set-A) line from its features
(no target is read; never a test line; the launcher needs the predicted line's lineage, so it must be in the split). Caches are redirected to a temporary directory (cold timing; the project's
cache/m2_inner, cache/m5_inner and cache/m5_outer are neither read nor written). Each model runs in its own
subprocess so peak RSS is its own. Output: results/timing/refit_full_dev.json.
"""
from __future__ import annotations

import os
import sys

if __name__ == "__main__":
    _n = sys.argv[sys.argv.index("--threads") + 1] if "--threads" in sys.argv else "4"
    for _v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[_v] = _n

import argparse  # noqa: E402
import datetime as dt  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "results" / "timing" / "refit_full_dev.json"
PREDICT_LINE = "ACH-000133"  # noncancerous DepMap line with every feature block (features only; no targets)


def one(name, threads):
    import pandas as pd
    import heldout_common as HC
    from essml.models import m2_local, m5_stack
    tmp = Path(tempfile.mkdtemp(prefix="essml_timing_"))
    m2_local.INNER_CACHE = tmp / "m2_inner"
    m5_stack.INNER_CACHE, m5_stack.OUTER_CACHE = tmp / "m5_inner", tmp / "m5_outer"
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(threads)
    except ImportError:
        pass
    sp = pd.read_csv(ROOT / "configs" / "splits.csv", index_col=0)
    if sp.loc[PREDICT_LINE, "role"] != "noncancerous":
        raise SystemExit(f"{PREDICT_LINE} must be a noncancerous (set-A) line, never a test line")
    cfg = json.loads((ROOT / "configs" / "runs" / f"{name}.json").read_text())
    f = HC.Finalist("timing", name, f"configs/runs/{name}.json", HC.config_hash(cfg), cfg)
    pop = HC.population(HC.Paths(), dry=False)
    t0 = time.perf_counter()
    r = HC.refit_predict(f, pop, HC.Paths(), predict_ids=[PREDICT_LINE],
                         params_extra={"cache": False} if cfg["model"] == "M5" else None)
    total = time.perf_counter() - t0
    rss, how = HC.peak_rss()
    rep = {"model": name, "n_train_ids": r["n_train_ids"], "threads": threads, "fit_s": round(r["fit_s"], 1),
           "predict_s": round(r["predict_s"], 1), "refit_s": round(r["fit_s"] + r["predict_s"], 1),
           "total_s_incl_target_load": round(total, 1), "peak_rss_mb": round(rss), "rss_method": how,
           "per_input_s": getattr(r["model"], "timing_", None), "predicted_line": PREDICT_LINE,
           "finite_prediction": bool(r["pred"].notna().all().all())}
    m3 = cfg.get("params", {}).get("inputs", [])
    rep["note"] = ("M5 average: fit stores the data; the three inputs are refitted inside predict (refit_s). M3 uses "
                   f"n_jobs={[i['params'].get('n_jobs') for i in m3 if i['model'] == 'M3']} worker processes (config)."
                   if cfg["model"] == "M5" else "fit_s includes the model's internal tuning (LOO lambda)")
    print(json.dumps(rep))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--models", default="m5_avg,m1_E_resB0")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--one", default=None, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.one:
        return one(a.one, a.threads)
    res = []
    for name in a.models.split(","):
        p = subprocess.run([sys.executable, __file__, "--one", name, "--threads", str(a.threads)],
                           capture_output=True, text=True, cwd=ROOT)
        if p.returncode != 0:
            raise SystemExit(f"{name} failed:\n{p.stderr[-3000:]}")
        res.append(json.loads(p.stdout.strip().splitlines()[-1]))
        print(res[-1])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"timestamp": dt.datetime.now().isoformat(timespec="seconds"),
                               "machine_note": "12 cores, 16 GB; nothing else running (checked before the run)",
                               "results": res}, indent=1))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
