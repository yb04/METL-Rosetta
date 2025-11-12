import pandas as pd
from collections import Counter
import yaml
from IPython.display import display

with open("config.yaml", 'r') as y:
    config = yaml.safe_load(y)
fn = config["dataPath"]


# 1) Schema-Sample nur für Spaltenliste und grundlegende Dtypes (keine parse_dates)
sample = pd.read_csv(fn, sep="\t", nrows=50_000, low_memory=False)


# 1) Header zuverlässig ermitteln (liest nur die Kopfzeile)
header_df = pd.read_csv(fn, sep="\t", nrows=0)
cols = header_df.columns.tolist()

# Die Datei beginnt mit einer leeren ersten Spalte (Index), die entfernen
if cols and cols[0] == "":
    cols = cols[1:]

# 2) Datums-/Zeitspalten komplett ausschließen
exclude_cols = [c for c in cols if any(k in c.lower() for k in ["date","time","timestamp"])]
use_cols = [c for c in cols if c not in exclude_cols]

# 3) Dtypes für Textspalten setzen (keine parse_dates!)
# mutations/pdb_fn/job_uuid sind Strings; Rest wird pro Chunk numerisch gecastet
string_like = [c for c in ["pdb_fn","mutations","job_uuid"] if c in use_cols]
dtype_map = {c: "string" for c in string_like}

rows = 0
mut_counter = Counter()
best_total_parts = []
pdb_mut_counts = Counter()

reader = pd.read_csv(
    fn,
    sep="\t",
    chunksize=500_000,
    usecols=use_cols,      # nur sauberen Satz an Spalten laden
    dtype=dtype_map,       # string-Dtype für die Textspalten
    low_memory=False
)

for chunk in reader:
    rows += len(chunk)

    # numerische Spalten konvertieren (alles außer string_like)
    num_cols = [c for c in chunk.columns if c not in string_like]
    chunk[num_cols] = chunk[num_cols].apply(pd.to_numeric, errors="coerce")

    # Mutationszählungen
    if "mutations" in chunk.columns:
        mut_series = chunk["mutations"].fillna("").astype("string")
        for ms in mut_series:
            for m in [x.strip() for x in ms.split(",") if x.strip()]:
                mut_counter[m] += 1

    # einzigartige Mutationslisten pro PDB
    if {"pdb_fn","mutations"}.issubset(chunk.columns):
        pdb_mut_counts.update(chunk.groupby("pdb_fn")["mutations"].nunique().to_dict())

    # Top-Designs nach total_score
    if "total_score" in chunk.columns:
        part = chunk.loc[:, ["pdb_fn","mutations","total_score"]].dropna()
        part = part.nsmallest(50, "total_score")
        best_total_parts.append(part)

# Ergebnisse
report_best_total = (
    pd.concat(best_total_parts).nsmallest(50, "total_score")
    if best_total_parts else pd.DataFrame()
)

top_mut = pd.Series(mut_counter).sort_values(ascending=False).head(30)
top_pdb = pd.Series(pdb_mut_counts).sort_values(ascending=False).head(30)

display({"rows_seen": rows})
display(top_mut.to_frame("count"))
display(top_pdb.to_frame("unique_mutation_sets"))
display(report_best_total)
