import gzip, bisect, collections

GFF = "data/annotation/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz"
chr_len, genes = {}, collections.defaultdict(list)
with gzip.open(GFF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        if p[2] == "chromosome":
            chr_len[p[0]] = int(p[4])
        elif p[2] == "gene":
            genes[p[0]].append((int(p[3]), int(p[4])))

def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

merged = {c: merge(v) for c, v in genes.items()}
snps = collections.defaultdict(list)
with open("data/work/snps_maf.bim") as f:
    for line in f:
        c, _, _, pos = line.split()[:4]
        snps[c].append(int(pos))

print(f"{'chr':>3} {'SNPs':>7} {'in genes %':>10} {'genic %':>8} {'ratio':>6} {'max SNP pos':>12} {'chr length':>11} {'within':>6}")
tot_n = tot_in = tot_genic = tot_len = 0
for c in [str(i) for i in range(1, 11)]:
    iv = merged.get(c, [])
    starts = [s for s, e in iv]
    n_in = 0
    for pos in snps[c]:
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and pos <= iv[i][1]:
            n_in += 1
    n = len(snps[c]); genic = sum(e - s + 1 for s, e in iv); L = chr_len[c]; mx = max(snps[c])
    print(f"{c:>3} {n:>7} {100*n_in/n:>10.1f} {100*genic/L:>8.1f} {(n_in/n)/(genic/L):>6.1f} {mx:>12} {L:>11} {str(mx <= L):>6}")
    tot_n += n; tot_in += n_in; tot_genic += genic; tot_len += L
print(f"all {tot_n:>7} {100*tot_in/tot_n:>10.1f} {100*tot_genic/tot_len:>8.1f} {(tot_in/tot_n)/(tot_genic/tot_len):>6.1f}")
