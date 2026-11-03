# Data description

## Source
Dryad dataset doi:10.5061/dryad.3j9kd51w9 (public domain): "Integrating SNP, RNASeq, and Phenotypic data from sorghum association and diversity panels for comprehensive genomic analysis" (Subash, Mangal, Shrestha, Xu, Schnable, Yang; University of Nebraska-Lincoln). Raw files are not stored in this repository (`data/` is git-ignored).

## Genotypes
File: `B.SDP_SAP_815_maf001_het10_imputed.recode.vcf.gz`
- 815 lines (Sorghum Association Panel plus Sorghum Diversity Panel), SNPs called from RNA-seq of mature leaf tissue against the BTx623 v5 reference (as stated in the dataset README; not independently verified).
- 346,222 variants: 277,724 SNPs (all biallelic) and 68,498 indels; chromosomes Chr01 to Chr10 only.
- Variants with MAF < 0.01 or > 10% heterozygous calls were removed by the data authors. The VCF header shows Beagle processing; no marker carries an imputation flag or INFO values, so imputed genotypes cannot be identified or filtered by imputation quality.
- No missing genotypes (0 of 282,170,930 calls); 1.07% of calls are heterozygous.
- Markers come from expressed regions, so genome coverage is uneven.

## Phenotypes
File: 2021 Lincoln, Nebraska field experiment (`2.SDP_SAP_f2021_field_phenotype`), high-nitrogen (HN) treatment only.
- Traits: plant height and days to flower, averaged per line over replicates (`metadata/pheno_2021_HN.tsv`, produced by `scripts/01_make_phenotype_table.py`).
- Plant height units: the dataset README says centimetres, but the values (median 0.94) indicate metres. Units are recorded here as metres (inferred).
- Two implausible plant height records (6.0 and 9.67) are treated as missing (rule: height > 5 m).
- Other years were not used: files differ in design and scale (for example the 2023 median height is 251 against 95.5 in 2020, and 2023 contains implausible values).

## Matching genotypes to phenotypes
Phenotype IDs contain spaces and varying case (for example "PI 152651"); the VCF uses names without spaces ("PI152651"). IDs are matched after removing spaces, hyphens and underscores and ignoring case, using the PI number, then the conversion ID, then the line name (scripts 00b and 01). In the 2021 file all 815 VCF lines have a matching entry.

## Known issues and open questions
- Replicate structure (Treatment x Rep) needs inspecting before modelling; tall lines appear mainly in HN replicates 4 and 5.
- One line (PI564163) is replicated in roughly 80 to 100 plots, so it is probably a repeated check; it is kept as a single line.
- Eight VCF sample names are not accession IDs (Border, Day, ElMota, RTAM2566, FeteritaGishesh, SanChiSan, Comb7078, Segaolane); Border is probably a field label and will be removed.
- The reference build (v5) and a matching gene annotation must be confirmed before candidate gene lookup (Week 15).
- The panel mixes SAP and SDP lines, so population structure is expected to be strong.

## Phenotype table summary
- Lines with plant height: 813; with days to flower: 813; with both: 813 (of 815 genotyped lines).
