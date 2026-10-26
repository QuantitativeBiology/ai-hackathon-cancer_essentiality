"""Download open external resources listed in configs/external_sources.json.

Each resource goes to data/raw_external/<name>/ (git-ignored). A manifest with URL, date,
size and sha256 is written to configs/external_manifest.json for DATA.md.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import subprocess


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", default="configs/external_sources.json")
    ap.add_argument("--out", default="data/raw_external")
    ap.add_argument("--manifest", default="configs/external_manifest.json")
    ap.add_argument("--only", nargs="*", help="subset of resource names")
    a = ap.parse_args()
    sources = {k: v for k, v in json.load(open(a.sources)).items() if not k.startswith("_")}
    man_path = pathlib.Path(a.manifest)
    manifest = json.load(open(man_path)) if man_path.exists() else {}
    for name, s in sources.items():
        if a.only and name not in a.only:
            continue
        d = pathlib.Path(a.out) / name
        d.mkdir(parents=True, exist_ok=True)
        fn = d / s["url"].rstrip("/").split("/")[-1]
        try:
            if not fn.exists():
                # curl: some providers return 403 to urllib's default client
                subprocess.run(["curl", "-sSfL", "-m", "1800", "-o", str(fn), s["url"]], check=True)
            manifest[name] = {"url": s["url"], "file": str(fn).replace("\\", "/"), "bytes": fn.stat().st_size,
                              "sha256": sha256(fn), "version": s["version"], "licence": s["licence"],
                              "retrieved": manifest.get(name, {}).get("retrieved", dt.date.today().isoformat())}
            print(f"OK   {name}: {fn.stat().st_size} bytes")
        except Exception as e:  # report and continue; failures are listed in the manifest
            fn.unlink(missing_ok=True)  # never leave a partial file that a later run would accept
            manifest[name] = {"url": s["url"], "error": repr(e)}
            print(f"FAIL {name}: {e!r}")
    json.dump(manifest, open(man_path, "w"), indent=1)


if __name__ == "__main__":
    main()
