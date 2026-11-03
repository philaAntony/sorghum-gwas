import csv, glob, gzip, re, collections

VCF = "data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz"
with gzip.open(VCF, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break

def norm(x):
    return re.sub(r"[\s_\-]", "", x).upper()

vcf_map = {norm(s): s for s in samples}

def name_type(s):
    u = s.upper()
    if u.startswith("PI"): return "PI"
    if u.startswith("SC"): return "SC"
    if u[:3] in ("BTX", "RTX"): return "BTx_RTx"
    return "other"

path = glob.glob("data/raw/2.*/*.csv")[0]
rows = list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")))
sap, conv = set(), set()
for r in rows:
    taxa = None
    for c in ("PINumber", "SorghumConversionID", "Name"):
        v = (r.get(c) or "").strip()
        if v and norm(v) in vcf_map:
            taxa = vcf_map[norm(v)]
            break
    if taxa is None:
        continue
    if (r.get("SAPID") or "").strip():
        sap.add(taxa)
    if (r.get("SorghumConversionID") or "").strip():
        conv.add(taxa)

with open("metadata/line_groups.tsv", "w") as out:
    out.write("taxa\tpanel\tname_type\n")
    for s in samples:
        out.write(f"{s}\t{'SAP' if s in sap else 'non-SAP'}\t{name_type(s)}\n")

keep = set(open("metadata/keep_samples.txt").read().split())
k = [s for s in samples if s in keep]
print("kept lines:", len(k))
print("panel x name type:", dict(collections.Counter(("SAP" if s in sap else "non-SAP", name_type(s)) for s in k)))
print("lines with a conversion ID:", sum(1 for s in k if s in conv))
