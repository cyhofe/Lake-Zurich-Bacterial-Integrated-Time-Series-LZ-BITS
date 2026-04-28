#!/usr/bin/env python3
"""
Folder-local KEGG BRITE mapper (no pandas, no absolute paths).

Looks in the SAME FOLDER as this script for:
- KeggAnnotation.txt
- Brite_Roots_Config.tsv.txt   (TSV with columns: major_group, brite_id, brite_name)

Writes to the SAME FOLDER:
- kegg_brite_master_long.tsv

Caches BRITE JSON to:
- brite_json/<brite_id>.json

Output columns:
KO | major_group | brite_id | brite_name | path | depth
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple
from urllib.request import urlopen, Request

# ----------------------------
# Settings
# ----------------------------
SLEEP_SECONDS = 0.40  # be polite to KEGG
KEGG_GET = "https://rest.kegg.jp/get"
KO_RE_EXACT = re.compile(r"^(K\d{5})$")
KO_RE_ANY = re.compile(r"\b(K\d{5})\b")
DELIM = " > "

# ----------------------------
# Locate files relative to script (not environment)
# ----------------------------
HERE = Path(__file__).resolve().parent

ANNOT_FILE = HERE / "KeggAnnotation.txt"
CONFIG_FILE = HERE / "Brite_Roots_Config.tsv.txt"
CACHE_DIR = HERE / "brite_json"
OUT_FILE = HERE / "kegg_brite_master_long.tsv"

CACHE_DIR.mkdir(exist_ok=True)


def die(msg: str) -> None:
    raise SystemExit(f"[ERROR] {msg}")


def extract_kos_from_annotation() -> Set[str]:
    if not ANNOT_FILE.exists():
        die(f"Missing input file: {ANNOT_FILE.name}")

    kos: Set[str] = set()
    with ANNOT_FILE.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            die(f"{ANNOT_FILE.name} must be a TSV with at least 2 columns")

        for row in reader:
            if len(row) < 2:
                continue
            val = row[1].strip()
            if not val or val == "-:-":
                continue
            ko = val.split(":", 1)[0].strip()
            if KO_RE_EXACT.match(ko):
                kos.add(ko)

    print(f"✓ Extracted {len(kos):,} unique KOs from {ANNOT_FILE.name}")
    return kos


def read_brite_config() -> List[Tuple[str, str, str]]:
    if not CONFIG_FILE.exists():
        die(f"Missing config file: {CONFIG_FILE.name}")

    rows: List[Tuple[str, str, str]] = []
    with CONFIG_FILE.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"major_group", "brite_id", "brite_name"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            die(
                f"{CONFIG_FILE.name} must be TSV with header: major_group, brite_id, brite_name\n"
                f"Found header: {reader.fieldnames}"
            )

        for r in reader:
            major = (r.get("major_group") or "").strip()
            brite_id = (r.get("brite_id") or "").strip()
            brite_name = (r.get("brite_name") or "").strip()
            if not brite_id:
                continue
            rows.append((major, brite_id, brite_name))

    print(f"✓ Loaded {len(rows):,} BRITE roots from {CONFIG_FILE.name}")
    return rows


def http_get_json(url: str) -> Dict:
    # KEGG sometimes blocks requests without a UA
    req = Request(url, headers={"User-Agent": "Mozilla/5.0 (KEGG-BRITE-mapper)"})
    with urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_brite_json(brite_id: str) -> Dict:
    cache = CACHE_DIR / f"{brite_id}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))

    url = f"{KEGG_GET}/br:{brite_id}/json"
    data = http_get_json(url)
    cache.write_text(json.dumps(data, indent=2), encoding="utf-8")
    time.sleep(SLEEP_SECONDS)
    return data


def walk_brite_tree(
    node: Dict,
    path_so_far: List[str],
    out_rows: List[Dict],
    major_group: str,
    brite_id: str,
    brite_name: str,
):
    name = (node.get("name") or "").strip()
    new_path = path_so_far + ([name] if name else [])

    m = KO_RE_ANY.search(name)
    if m:
        out_rows.append(
            {
                "KO": m.group(1),
                "major_group": major_group,
                "brite_id": brite_id,
                "brite_name": brite_name,
                "path": DELIM.join(new_path),
                "depth": str(len(new_path)),
            }
        )

    for child in node.get("children") or []:
        walk_brite_tree(child, new_path, out_rows, major_group, brite_id, brite_name)


def write_tsv(rows: List[Dict]) -> None:
    with OUT_FILE.open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["KO", "major_group", "brite_id", "brite_name", "path", "depth"]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    my_kos = extract_kos_from_annotation()
    cfg = read_brite_config()

    all_rows: List[Dict] = []

    for major_group, brite_id, brite_name in cfg:
        print(f"Processing {brite_id} ({brite_name})")
        data = fetch_brite_json(brite_id)
        walk_brite_tree(data, [], all_rows, major_group, brite_id, brite_name)

    # Filter to observed KOs
    filtered = [r for r in all_rows if r["KO"] in my_kos]

    # Deduplicate exact duplicates
    seen = set()
    uniq: List[Dict] = []
    for r in filtered:
        key = (r["KO"], r["major_group"], r["brite_id"], r["brite_name"], r["path"], r["depth"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)

    write_tsv(uniq)

    # Simple summary
    kos_mapped = {r["KO"] for r in uniq}
    from collections import Counter
    c = Counter([r["KO"] for r in uniq])
    multi = sum(1 for _, v in c.items() if v > 1)

    print("\n✅ DONE")
    print(f"Output: {OUT_FILE.name}")
    print(f"Rows: {len(uniq):,}")
    print(f"Unique KOs mapped: {len(kos_mapped):,}")
    print(f"Multi-mapped KOs: {multi:,}")
    print(f"Cache: {CACHE_DIR.name}/")


if __name__ == "__main__":
    main()
