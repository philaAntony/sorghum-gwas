import csv, gzip, collections

GFF = "data/annotation/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz"
genes = collections.defaultdict(list)
with gzip.open(GFF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        genes[p[0]].append({"id": a["ID"].replace(".v5.1", ""), "start": int(p[3]),
                            "end": int(p[4]), "strand": p[6], "biotype": a.get("biotype", "")})

regions = list(csv.DictReader(open("results/candidate_regions.tsv"), delimiter="\t"))
rows, bad_id = [], 0
for r in regions:
    trait, chrom, lead = r["trait"], r["chr"], r["lead_snp"]
    lead_pos = int(lead.split(":")[1])
    rs, re_ = int(r["region_start"]), int(r["region_end"])
    ss, se = int(r["sig_start"]), int(r["sig_end"])
    snps = []
    with gzip.open(f"results/gwas/{trait}_MLM_K_PC3.tsv.gz", "rt") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d["chr"] == chrom and rs - 200000 <= int(d["pos"]) <= re_ + 200000:
                snps.append((int(d["pos"]), float(d["p"])))
    for g in genes[chrom]:
        if g["end"] < rs or g["start"] > re_:
            continue
        if int(g["id"][6:9]) != int(chrom):
            bad_id += 1
        inside = [p for pos, p in snps if g["start"] <= pos <= g["end"]]
        dist = 0 if g["start"] <= lead_pos <= g["end"] else min(abs(g["start"] - lead_pos), abs(g["end"] - lead_pos))
        rows.append({"trait": trait, "chr": chrom, "gene": g["id"], "start": g["start"], "end": g["end"],
                     "strand": g["strand"], "biotype": g["biotype"], "dist_to_lead_bp": dist,
                     "in_sig_range": int(g["end"] >= ss and g["start"] <= se),
                     "n_snps": len(inside), "best_p": min(inside) if inside else ""})

rows.sort(key=lambda x: (x["trait"], int(x["chr"]), x["start"]))
cols = list(rows[0].keys())
with open("results/candidate_genes.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

print("gene IDs whose chromosome number disagrees with the region:", bad_id)
by = collections.defaultdict(list)
for x in rows:
    by[(x["trait"], x["chr"])].append(x)
for (t, c), v in by.items():
    withsnp = [x for x in v if x["n_snps"]]
    print(f"\n== {t}, chr{c}: {len(v)} genes in region, {len(withsnp)} with SNPs inside, "
          f"{sum(x['dist_to_lead_bp'] == 0 for x in v)} containing the lead SNP")
    for x in sorted(withsnp, key=lambda x: x["best_p"])[:8]:
        print(f"  {x['gene']}  {x['start']}-{x['end']} ({x['strand']})  SNPs {x['n_snps']}  best p {x['best_p']:.2e}  distance to lead {x['dist_to_lead_bp']}")
