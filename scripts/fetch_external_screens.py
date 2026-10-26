"""Re-download external CRISPR screens and expression sources listed in configs/external_screens.json.

Every file goes to data/<root>/<deposit>/<name> (git-ignored), one new directory per deposit.
URLs come only from the config (direct `url` or a `template` + fields). Downloads use curl
(`-sSfL`, resume with `-C -`, retries) because some servers reject python urllib. Nothing
downloaded is executed or imported. A manifest with URL, bytes, sha256 and retrieval date is
written to configs/external_screens_manifest.json; an `md5` in the config is verified.

Usage:  python scripts/fetch_external_screens.py [--only DEPOSIT ...] [--verify-only] [--force]
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stub(acc):
    """GEO directory stub: GSE128210 -> GSE128nnn, GSE1234 -> GSEnnn."""
    prefix, num = acc[:3], acc[3:]
    return prefix + (num[:-3] if len(num) > 3 else "") + "nnn"


def resolve_url(f, templates):
    if "url" in f:
        return f["url"]
    fields = dict(f)
    if "gse" in f:
        fields["gse_stub"] = stub(f["gse"])
    if "gsm" in f:
        fields["gsm_stub"] = stub(f["gsm"])
    return templates[f["template"]].format(**fields)


def curl(url, dest, retries=6, timeout=3600):
    """Download with resume; retry on transient failures. Returns None or raises."""
    last = None
    for attempt in range(retries):
        r = subprocess.run(["curl", "-sSfL", "-C", "-", "-m", str(timeout), "-o", str(dest), url],
                           capture_output=True, text=True)
        if r.returncode == 0:
            return
        last = f"curl exit {r.returncode}: {r.stderr.strip()}"
        if r.returncode == 33:  # server refuses ranges: restart from scratch
            dest.unlink(missing_ok=True)
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(last)


def looks_like_html(path):
    with open(path, "rb") as fh:
        head = fh.read(512).lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "configs/external_screens.json"))
    ap.add_argument("--manifest", default=str(ROOT / "configs/external_screens_manifest.json"))
    ap.add_argument("--only", nargs="*", help="subset of deposit ids")
    ap.add_argument("--verify-only", action="store_true", help="hash existing files, download nothing")
    ap.add_argument("--force", action="store_true", help="re-download even if the file exists")
    a = ap.parse_args()

    cfg = json.load(open(a.config, encoding="utf-8"))
    templates = cfg["url_templates"]
    man_path = pathlib.Path(a.manifest)
    manifest = json.load(open(man_path, encoding="utf-8")) if man_path.exists() else {}
    n_ok = n_fail = 0
    for dep_id, dep in cfg["deposits"].items():
        if a.only and dep_id not in a.only:
            continue
        d = ROOT / "data" / dep["root"] / dep_id
        d.mkdir(parents=True, exist_ok=True)
        for f in dep["files"]:
            url = resolve_url(f, templates)
            dest = d / f["name"]
            key = f"{dep_id}/{f['name']}"
            try:
                if a.force and dest.exists() and not a.verify_only:
                    dest.unlink()
                if not dest.exists():
                    if a.verify_only:
                        raise FileNotFoundError(dest)
                    curl(url, dest)
                if looks_like_html(dest) and not f["name"].endswith((".html", ".htm")):
                    raise RuntimeError("server returned an HTML page (bot check or error page), not the file")
                if "md5" in f and digest(dest, "md5") != f["md5"]:
                    raise RuntimeError("md5 mismatch")
                prev = manifest.get(key, {})
                manifest[key] = {"url": url, "file": str(dest.relative_to(ROOT)).replace("\\", "/"),
                                 "bytes": dest.stat().st_size, "sha256": digest(dest, "sha256"),
                                 "licence": dep.get("licence"), "source": dep.get("source"),
                                 "retrieved": prev.get("retrieved", dt.date.today().isoformat())}
                n_ok += 1
                print(f"OK   {key} {dest.stat().st_size}")
            except Exception as e:  # report, keep going; never keep a partial/invalid file
                if not a.verify_only and dest.exists() and not isinstance(e, FileNotFoundError):
                    dest.unlink()
                manifest[key] = {"url": url, "error": repr(e)}
                n_fail += 1
                print(f"FAIL {key}: {e!r}")
    man_path.write_text(json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    print(f"{n_ok} ok, {n_fail} failed -> {man_path}")


if __name__ == "__main__":
    main()
