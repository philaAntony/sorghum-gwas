# GWAS theory notes 

1. Linkage disequilibrium (LD) and causal variants
LD is the non-random association between alleles at different loci. Variants close together on a chromosome tend to be inherited together because recombination has had little chance to separate them. It is usually measured by r² (0 = independent, 1 = perfectly correlated).
A GWAS tests each SNP individually. If a causal variant sits in a haplotype block with 50 other SNPs, all of them will correlate with the trait because they correlate with the causal variant. The result is a peak of significant SNPs, and the most significant one is often just the best “tag,” not the cause. Genotyping noise, allele frequency differences, and sample size all affect which one comes out on top. Fine-mapping, functional annotation, and experimental validation are needed to pin down the real causal variant. My markers come from RNA-seq, so the causal variant may not be among them, and LD decay decides how wide my candidate window is.
2. Population structure
Population structure is systematic differences in allele frequencies between subgroups of your sample, due to different ancestry, demographic history, or geography.
False association example: Suppose your sample mixes two ancestral groups. Group A has a higher average height (for genetic or environmental reasons), and also happens to have a higher frequency of allele T at a SNP unrelated to height. Pooled together, T carriers are taller on average, and you get a significant association, even though T has no effect within either group. The classic textbook example is the chopsticks-use “gene”: any allele more common in East Asian ancestry would correlate with chopstick use.
The panel combines the SAP (adapted lines with introgressed dwarfing genes) and the diversity panel, so height and allele frequencies can differ between groups.
3. Kinship matrix
A kinship matrix (K) is an n × n matrix, where n is the number of individuals, giving the pairwise genetic relatedness between every pair of samples. It captures both close relatives and subtler ancestry sharing.
Estimating from SNPs (a common approach, VanRaden’s method):
1.	Build a genotype matrix M (individuals × SNPs) coded 0/1/2 copies of the alternate allele.
2.	Center each SNP by its mean (2p, where p is the allele frequency), giving Z.
3.	Compute K = ZZᵀ / (2 Σ pᵢ(1 − pᵢ)), where the sum is over SNPs.
Diagonal entries reflect inbreeding or self-relatedness; off-diagonals measure sharing between individuals. Many pipelines use leave-one-chromosome-out (LOCO) kinship to avoid the tested SNP’s region inflating its own signal.
4. Fixed and random effects in a mixed model
The model is typically y = Xβ + Zu + ε:
•	Fixed effects (β): Effects you estimate directly and treat as specific, known quantities: the SNP being tested, and often top principal components. These are the things you want an effect size and p-value for.
•	Random effects (u): The polygenic background, i.e., the combined small effects of all genetic variation across the genome. It is modeled as u ~ N(0, σ²ₘK). This soaks up structure and relatedness.
•	Residual (ε): Independent noise, ε ~ N(0, σ²ₑI).
The ratio of σ²ₘ to total variance gives a heritability estimate. Testing the SNP as a fixed effect while the random effect absorbs background relatedness is what protects against false positives.
5. QQ plot
A QQ plot compares the observed −log₁₀(p) values (y-axis, sorted) against the values expected under the null (x-axis, from a uniform distribution). Under the null, points follow the diagonal y = x.
•	Well-calibrated: Points hug the diagonal, then lift off only at the extreme tail, reflecting true associations.
•	Uncorrected structure (inflation): The entire curve shifts above the diagonal, starting early, because many SNPs show spurious associations. This is quantified by the genomic inflation factor λ, the median observed χ² divided by the expected median (0.456). λ ≈ 1 is good; λ > ~1.05–1.1 suggests uncorrected structure or relatedness.
Caveat: for highly polygenic traits, genuine signal can also cause mild inflation, so λ alone doesn’t prove structure is the problem. LD score regression can help distinguish them.
6. Significance threshold
Testing ~10⁶ SNPs at α = 0.05 would give ~50,000 false positives, so you must correct for multiple testing.
•	Bonferroni: α / m. It is derived from the effective number of independent common-variant tests in the genome, since LD means SNPs aren’t truly independent.
•	Effective number of tests (Meff): Estimate the number of independent tests accounting for LD, then apply Bonferroni. 
•	Permutation: Shuffle phenotypes many times and record the minimum p-value each time. The 5th percentile of those minima is your empirical threshold. It is accurate but computationally heavy.
•	FDR (Benjamini-Hochberg): Controls the expected proportion of false discoveries rather than the chance of any false positive. It is more powerful, and common in studies like eQTL mapping, but less often used for headline GWAS claims.
•	Suggestive threshold: Often ~1 × 10⁻⁵ to flag loci worth following up even if they don’t reach genome-wide significance. Studies with many rare variants or sequencing data may use stricter thresholds (e.g., ~5 × 10⁻⁹) because there are more independent tests.

