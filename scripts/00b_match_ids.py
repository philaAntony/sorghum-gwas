import csv, glob, gzip, re, collections

VCF = "data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz"
with gzip.open(VCF, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break

def norm(x):
    return re.sub(r"[\s_\-]", "", x).upper()

def kind(s):
    u = s.upper()
    if u.startswith("PI"): return "PI"
    if u.startswith("SC"): return "SC"
    if u[:3] in ("BTX", "RTX"): return "BTx/RTx"
    return "other"

vcf_map = collections.defaultdict(list)
for s in samples:
    vcf_map[norm(s)].append(s)
print("VCF samples:", len(samples), "| unique after normalising:", len(vcf_map))
print("name types:", dict(collections.Counter(kind(s) for s in samples)))
print("examples of 'other':", [s for s in samples if kind(s) == "other"][:12])

def read(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

datasets = [
    ("2020 Lincoln", glob.glob("data/raw/1.*/*.csv")[0], ["SorghumAccession", "SorghumName"], "Treatment", ["PlantHeight", "DaysToBloom"]),
    ("2021 Lincoln", glob.glob("data/raw/2.*/*.csv")[0], ["PINumber", "SorghumConversionID", "Name"], "Treatment", ["Plant Height", "DaysToFlower"]),
    ("2022 Alabama", glob.glob("data/raw/3.*/*.csv")[0], ["Genotype"], "Treatment", ["HeadingDate"]),
    ("2023 Lincoln", "data/raw/4A.SAP_2023_field_phenotype.csv.gz", ["geno"], "Treatment_Type", ["PlantHeight", "DaysToBloom"]),
    ("2023 greenhouse", "data/raw/4B.SAP_2023_GH_phenotype.csv.gz", ["Taxa"], None, ["PH.HN", "PH.LN"]),
    ("2024 Lincoln", "data/raw/5.SAP_2024_field_Phenotype.csv.gz", ["geno"], "experiment", ["Plant Height 1", "Plant Height 2"]),
]

for name, path, id_cols, tcol, traits in datasets:
    rows = read(path)
    matched, n_match = set(), 0
    unmatched = collections.Counter()
    counts = collections.Counter()
    for r in rows:
        hit = None
        for c in id_cols:
            v = (r.get(c) or "").strip()
            if v and v.upper() != "NA" and norm(v) in vcf_map:
                hit = norm(v)
                break
        if hit:
            n_match += 1
            matched.add(hit)
            tr = (r.get(tcol) or "?") if tcol else "all"
            for t in traits:
                if num(r.get(t)) is not None:
                    counts[(tr, t)] += 1
        else:
            first = next(((r.get(c) or "").strip() for c in id_cols if (r.get(c) or "").strip()), "")
            unmatched[first] += 1
    print(f"\n== {name}: {len(rows)} rows, {n_match} matched, {len(matched)} unique lines matched")
    for k in sorted(counts):
        print("   non-missing", k[1], "|", k[0], ":", counts[k])
    print("   most common unmatched IDs:", unmatched.most_common(8))
