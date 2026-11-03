import csv, glob, gzip

VCF = "data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz"
with gzip.open(VCF, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            samples = set(line.rstrip("\n").split("\t")[9:])
            break
print("VCF samples:", len(samples))

def check(path, cols):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    print("\n" + path.split("/")[-1], "| rows:", len(rows))
    for c in cols:
        vals = {r[c].strip() for r in rows if r.get(c) and r[c].strip()}
        hit = vals & samples
        print(f"  {c}: {len(vals)} unique, {len(hit)} found in VCF, e.g. {sorted(vals)[:3]}")

check(glob.glob("data/raw/1.*/*.csv")[0], ["SorghumAccession", "SorghumName", "SNPDataID"])
check(glob.glob("data/raw/2.*/*.csv")[0], ["PINumber", "SAPID", "SorghumConversionID", "Name"])
check(glob.glob("data/raw/3.*/*.csv")[0], ["Genotype"])
check("data/raw/4A.SAP_2023_field_phenotype.csv.gz", ["geno"])
check("data/raw/4B.SAP_2023_GH_phenotype.csv.gz", ["Taxa"])
check("data/raw/5.SAP_2024_field_Phenotype.csv.gz", ["geno"])
