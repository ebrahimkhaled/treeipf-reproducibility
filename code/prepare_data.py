# Build the two analysis files used by 02_run_clinical_studies.R from the raw downloads.
#
#   python code/prepare_data.py [RAW_DIR]
#
# RAW_DIR must contain
#   ist/IST_corrected.csv                     International Stroke Trial, corrected v2 (Edinburgh DataShare)
#   physionet2012/set-{a,b,c}.tar.gz          PhysioNet/CinC Challenge 2012 v1.0.0 records
#   physionet2012/Outcomes-{a,b,c}.txt        and their outcome files
# and optionally MANIFEST_SHA256.txt, against which the raw files are checked.
#
# Writes data/ist_dd6.parquet, data/physionet_ihd.parquet and data/data_provenance.csv.
# The feature definitions and exclusions follow the builder these files were first made
# with (paper_tfm_calibration/rebuild/code/realdata.py); only the columns this paper uses are kept.

import argparse
import csv
import hashlib
import io
import math
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_RAW = r"D:\Research\paper_tfm_calibration\rebuild\data\raw"
DATA = Path(__file__).resolve().parents[1] / "data"

N_COUNTRIES = 12


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_manifest(raw, files):
    man = raw / "MANIFEST_SHA256.txt"
    if not man.exists():
        print("  no MANIFEST_SHA256.txt in the raw folder; hashes not checked")
        return
    known = {}
    for line in man.read_text(encoding="utf8").splitlines():
        if line.strip() and not line.startswith("#"):
            digest, rel = line.split("  ", 1)
            known[rel.strip()] = digest
    for p in files:
        rel = p.relative_to(raw).as_posix()
        if rel not in known:
            print(f"  {rel}: not in the manifest")
        elif sha256(p) != known[rel]:
            sys.exit(f"{rel}: SHA256 does not match MANIFEST_SHA256.txt")
    print(f"  {len(files)} raw files match MANIFEST_SHA256.txt")


# ---------------------------------------------------------------------------
# International Stroke Trial
# ---------------------------------------------------------------------------

IST_OURS = ["age", "female", "rsbp", "rconsc", "ratrial", "ct", "rdef1", "rdef2"]


def yes_no(s):
    return s.astype(str).str.strip().str.upper().map({"Y": 1.0, "N": 0.0})


def upper(s):
    return s.astype(str).str.strip().str.upper()


def build_ist(raw, flow):
    path = raw / "ist" / "IST_corrected.csv"
    # the file is not utf-8; latin-1 reads every byte and the columns used here are plain ascii
    r = pd.read_csv(path, low_memory=False, encoding="latin-1")
    r.columns = [c.strip().upper() for c in r.columns]
    if len(r) != 19435:
        sys.exit(f"IST: expected 19435 patients, found {len(r)}")

    d = pd.DataFrame({"id": np.arange(len(r))})       # row position in IST_corrected.csv, from 0
    d["country"] = pd.to_numeric(r["CNTRYNUM"], errors="coerce")
    d["age"] = pd.to_numeric(r["AGE"], errors="coerce")
    d["female"] = upper(r["SEX"]).map({"F": 1.0, "M": 0.0})
    d["rsbp"] = pd.to_numeric(r["RSBP"], errors="coerce")
    d["rconsc"] = upper(r["RCONSC"]).where(lambda s: s.isin(["F", "D", "U"]))
    d["ratrial"] = yes_no(r["RATRIAL"])               # blank for the pilot-phase patients
    # ct: none = no scan before randomisation; infarct = scan showed a visible infarct
    rct, rvis = yes_no(r["RCT"]), yes_no(r["RVISINF"])
    ct = np.where(rct == 0, "none", np.where(rvis == 1, "infarct", np.where(rvis == 0, "no_infarct", None)))
    d["ct"] = pd.Series(ct, dtype=object).where(rct.notna(), None)
    # face and arm/hand deficit: Y = 1, N = 0, C (cannot assess, e.g. drowsy patient) = 0
    for i in (1, 2):
        d[f"rdef{i}"] = upper(r[f"RDEF{i}"]).map({"Y": 1.0, "N": 0.0, "C": 0.0})

    # outcome: OCCODE 1 dead, 2 dependent, 3 not recovered, 4 recovered; 0 and 9 unknown
    occ = pd.to_numeric(r["OCCODE"], errors="coerce")
    d["y"] = np.where(occ.isin([1, 2]), 1.0, np.where(occ.isin([3, 4]), 0.0, np.nan))

    def step(name, t, note=""):
        events = int(t["y"].sum()) if t["y"].notna().all() else None
        flow.append({"dataset": "IST", "step": name, "n": len(t), "events": events, "note": note})

    step("raw records in IST_corrected.csv", d)
    d = d.loc[d["y"].notna()].copy()
    d["y"] = d["y"].astype(int)
    step("6-month outcome known (OCCODE 1-4)", d, "OCCODE 0 or 9 excluded")
    miss = d[IST_OURS].isna().sum()
    d = d.loc[d[IST_OURS + ["country"]].notna().all(axis=1)].reset_index(drop=True)
    step("complete on the 8 predictors used here", d,
         "missing: " + ", ".join(f"{v} {n}" for v, n in miss.items() if n)
         + " (RATRIAL was not collected in the pilot phase, the first 984 patients)")

    # same rule as the R script: largest 12 countries in the complete-case table
    sizes = d.groupby("country").size().sort_values(ascending=False, kind="stable")
    if sizes.iloc[N_COUNTRIES - 1] == sizes.iloc[N_COUNTRIES]:
        sys.exit("IST: tie at the 12th country")
    top = d.loc[d["country"].isin(sizes.index[:N_COUNTRIES])]
    step(f"top {N_COUNTRIES} countries by size (analysed in the paper)", top,
         "country numbers " + " ".join(str(c) for c in sizes.index[:N_COUNTRIES]) + "; selected in 02_run_clinical_studies.R")

    # the file keeps every country; 02_run_clinical_studies.R picks the 12
    return d[["id", "country", "y"] + IST_OURS]


# ---------------------------------------------------------------------------
# PhysioNet/CinC Challenge 2012
# ---------------------------------------------------------------------------

DESCRIPTORS = {"RecordID", "Age", "Gender", "Height", "ICUType"}
# values outside these ranges are recording errors and treated as missing
RANGE = {"HR": (1, 300), "SysABP": (1, 300), "NISysABP": (1, 300), "Temp": (25, 45),
         "GCS": (3, 15), "BUN": (0, 300), "Creatinine": (0.05, 30), "WBC": (0.05, 200)}
FEATURES = {"hr_max": (("HR",), max), "sbp_min": (("SysABP", "NISysABP"), min),
            "temp_max": (("Temp",), max), "gcs_min": (("GCS",), min),
            "bun_max": (("BUN",), max), "creat_max": (("Creatinine",), max),
            "wbc_max": (("WBC",), max)}
PHY_OURS = ["age", "female", "hr_max", "sbp_min", "temp_max", "gcs_min",
            "bun_max", "creat_max", "wbc_max", "mechvent"]


def to_float(x):
    try:
        return float(x)
    except ValueError:
        return math.nan


def read_record(text):
    # each record holds the first 48 hours of one ICU stay: Time (hh:mm from admission), Parameter, Value
    desc, series = {}, {}
    rows = csv.reader(io.StringIO(text))
    next(rows)
    for time, prm, val in rows:
        v = to_float(val)
        if time == "00:00" and prm not in desc:
            desc[prm] = v
        if prm in DESCRIPTORS or (time == "00:00" and prm == "Weight"):
            continue
        lo, hi = RANGE.get(prm, (-math.inf, math.inf))
        if math.isnan(v) or not lo <= v <= hi:
            continue
        series.setdefault(prm, []).append(v)

    age, sex = desc.get("Age", -1), desc.get("Gender", -1)
    row = {"id": int(desc["RecordID"]),
           "age": age if age >= 0 else math.nan,
           "female": {0.0: 1.0, 1.0: 0.0}.get(sex, math.nan)}     # challenge coding: 0 female, 1 male
    for name, (prms, how) in FEATURES.items():
        vals = [v for p in prms for v in series.get(p, [])]
        row[name] = how(vals) if vals else math.nan
    # MechVent is only charted (as 1) while the patient is ventilated
    row["mechvent"] = float(any(v > 0 for v in series.get("MechVent", [])))
    return row


def build_physionet(raw, flow):
    d0 = raw / "physionet2012"
    recs, per_set = [], {}
    for s in "abc":
        with tarfile.open(d0 / f"set-{s}.tar.gz", "r:gz") as tf:
            members = [m for m in tf.getmembers() if m.isfile() and m.name.endswith(".txt")]
            per_set[s] = len(members)
            for m in members:
                recs.append(read_record(tf.extractfile(m).read().decode("utf8")))
    feats = pd.DataFrame(recs)
    oc = pd.concat([pd.read_csv(d0 / f"Outcomes-{s}.txt") for s in "abc"], ignore_index=True)
    if len(feats) != 12000 or len(oc) != 12000:
        sys.exit(f"PhysioNet: expected 12000 records, got {len(feats)} and {len(oc)}")

    d = feats.merge(oc[["RecordID", "In-hospital_death"]].rename(
        columns={"RecordID": "id", "In-hospital_death": "y"}), on="id", how="inner", validate="one_to_one")
    step = lambda name, t, note="": flow.append(
        {"dataset": "PhysioNet2012", "step": name, "n": len(t), "events": int(t["y"].sum()), "note": note})
    step("ICU stays in sets A, B, C", d, ", ".join(f"set {s.upper()} {n}" for s, n in per_set.items()))
    d = d.loc[d["y"].isin([0, 1])].sort_values("id").reset_index(drop=True)
    d["y"] = d["y"].astype(int)
    step("in-hospital outcome known", d)
    miss = d[PHY_OURS].isna().sum()
    cc = d.loc[d[PHY_OURS].notna().all(axis=1)]
    step("complete on the 10 predictors used here (analysed in the paper)", cc,
         "missing: " + ", ".join(f"{v} {n}" for v, n in miss.items() if n))

    # the file keeps all stays; 02_run_clinical_studies.R takes the complete cases
    return d[["id", "y"] + PHY_OURS]


def main():
    ap = argparse.ArgumentParser(description="Build the analysis files from the raw downloads.")
    ap.add_argument("raw_dir", nargs="?", default=DEFAULT_RAW)
    raw = Path(ap.parse_args().raw_dir)

    check_manifest(raw, [raw / "ist" / "IST_corrected.csv"]
                   + [raw / "physionet2012" / f"set-{s}.tar.gz" for s in "abc"]
                   + [raw / "physionet2012" / f"Outcomes-{s}.txt" for s in "abc"])
    flow = []
    ist = build_ist(raw, flow)
    phy = build_physionet(raw, flow)
    DATA.mkdir(exist_ok=True)
    ist.to_parquet(DATA / "ist_dd6.parquet", index=False)
    phy.to_parquet(DATA / "physionet_ihd.parquet", index=False)

    flow = pd.DataFrame(flow)
    flow["excluded"] = -flow.groupby("dataset")["n"].diff()
    flow["excluded"] = flow["excluded"].astype("Int64")
    flow["events"] = flow["events"].astype("Int64")
    flow = flow[["dataset", "step", "n", "excluded", "events", "note"]]
    flow.to_csv(DATA / "data_provenance.csv", index=False)
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(flow.drop(columns="note").to_string(index=False))
    print(f"wrote ist_dd6.parquet {ist.shape}, physionet_ihd.parquet {phy.shape}, data_provenance.csv")


if __name__ == "__main__":
    main()
