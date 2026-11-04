import csv, gzip, math, collections

GFF = "data/annotation/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz"
KNOWN = {"Dw1": "Sobic.009G229800", "Dw2": "Sobic.006G067700", "Dw3": "Sobic.007G163800",
         "Ma1": "Sobic.006G057866", "Ma6": "Sobic.006G004400"}
genes = {}
with gzip.open(GFF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        genes[a["ID"].replace(".v5.1", "")] = (p[0], int(p[3]), int(p[4]), p[6])

tabs = {}
for tr in ("log_plant_height", "days_to_flower"):
    d = collections.defaultdict(list)
    with gzip.open(f"results/gwas/{tr}_MLM_K_PC3.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                d[r["chr"]].append((int(r["pos"]), float(r["p"])))
            except ValueError:
                pass
    tabs[tr] = d

def best(tr, chrom, a, b):
    v = [p for pos, p in tabs[tr][chrom] if a <= pos <= b]
    return (len(v), -math.log10(max(min(v), 1e-300))) if v else (0, None)

regs = list(csv.DictReader(open("results/candidate_regions.tsv"), delimiter="\t"))
print("thresholds: Bonferroni -log10 6.52, 1/m -log10 5.22")
for name, gid in KNOWN.items():
    if gid not in genes:
        near = sorted(g for g in genes if g.startswith(gid[:-2]))
        print(f"\n{name} {gid}: NOT in v5.1 annotation; genes with the same ID prefix: {near}")
        continue
    chrom, s, e, st = genes[gid]
    print(f"\n{name} {gid}: chr{chrom}:{s}-{e} ({st})")
    for r in regs:
        if r["chr"] == chrom:
            a, b = int(r["region_start"]), int(r["region_end"])
            d = 0 if s <= b and e >= a else min(abs(s - b), abs(e - a))
            print(f"   region {r['trait']} chr{chrom}:{a}-{b}: distance {d} bp")
    for tr in tabs:
        cells = []
        for pad in (0, 100000, 500000):
            n, b = best(tr, chrom, s - pad, e + pad)
            cells.append("none" if b is None else f"{b:.2f} (n={n})")
        print(f"   {tr:17} best -log10 p in gene / +-100 kb / +-500 kb: {cells[0]} / {cells[1]} / {cells[2]}")
