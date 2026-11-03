import csv, glob, gzip, re, statistics, collections

VCF = "data/raw/B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz"
with gzip.open(VCF, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break

def norm(x):
    return re.sub(r"[\s_\-]", "", x).upper()

vcf_map = {norm(s): s for s in samples}

path = glob.glob("data/raw/2.*/*.csv")[0]
with open(path, newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

ph, dfl = collections.defaultdict(list), collections.defaultdict(list)
MAX_HEIGHT = 5.0  # metres; larger values treated as implausible
excluded = []
via, unmatched = collections.Counter(), collections.Counter()
for r in rows:
    if (r.get("Treatment") or "").strip() != "HN":
        continue
    taxa = None
    for c in ("PINumber", "SorghumConversionID", "Name"):
        v = (r.get(c) or "").strip()
        if v and norm(v) in vcf_map:
            taxa = vcf_map[norm(v)]
            via[c] += 1
            break
    if taxa is None:
        unmatched[(r.get("PINumber") or r.get("Name") or "").strip()] += 1
        continue
    a, b = num(r.get("Plant Height")), num(r.get("DaysToFlower"))
    if a is not None and a > MAX_HEIGHT:
        excluded.append((taxa, a))
        a = None
    if a is not None:
        ph[taxa].append(a)
    if b is not None:
        dfl[taxa].append(b)

def mean(v):
    return f"{statistics.mean(v):.2f}" if v else "NA"

with open("metadata/pheno_2021_HN.tsv", "w") as out:
    out.write("taxa\tplant_height\tn_ph\tdays_to_flower\tn_df\n")
    for s in samples:
        out.write(f"{s}\t{mean(ph[s])}\t{len(ph[s])}\t{mean(dfl[s])}\t{len(dfl[s])}\n")

print("HN rows matched via:", dict(via))
print("HN rows unmatched:", sum(unmatched.values()), "| top:", unmatched.most_common(6))
for name, d in (("plant height", ph), ("days to flower", dfl)):
    means = [statistics.mean(v) for v in d.values() if v]
    print(f"\n{name}: {len(means)} of {len(samples)} VCF lines have data")
    print("  replicates per line:", dict(sorted(collections.Counter(len(v) for v in d.values()).items())))
    print(f"  line means: mean {statistics.mean(means):.1f}, sd {statistics.stdev(means):.1f}, min {min(means):.1f}, max {max(means):.1f}")
both = [s for s in samples if ph[s] and dfl[s]]
print("\nlines with both traits:", len(both))
print("VCF lines with no plant height:", [s for s in samples if not ph[s]][:15])

print("heights excluded as implausible (> 5):", excluded)
