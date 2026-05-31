# BileActome Functional Gene Quantification Pipeline

This repository provides a reusable six-step workflow for annotating bile-acid-metabolizing genes in a non-redundant gene catalog using custom HMM profiles, extracting matched coding sequences, mapping metagenomic reads back to the matched functional CDSs, and generating gene-level and function-level abundance matrices.

The workflow is designed for BileActome-style HMM resources but can be reused for any functional HMM database where users want to quantify gene families from metagenomic read data.

---

## 1. Overview

### Purpose

The pipeline quantifies functional genes from a non-redundant gene catalog using a custom HMM database.

The main workflow is:

1. Search predicted proteins against custom HMM profiles using `hmmscan`.
2. Extract matched CDS sequences from the corresponding nucleotide gene catalog.
3. Map clean paired-end metagenomic reads to the matched functional CDSs using Bowtie2.
4. Calculate gene-level RPKM abundance for each sample.
5. Merge all sample-level RPKM tables and add functional annotation.
6. Aggregate gene-level abundance to function-level abundance and generate relative abundance and CLR matrices.

### Input concept

The pipeline assumes that the user has:

- A non-redundant protein catalog.
- A matching non-redundant CDS catalog.
- A directory of HMM profiles representing target functional gene families.
- Clean paired-end metagenomic reads for each sample.

### Output concept

The final outputs include:

- A gene-to-function mapping table.
- A FASTA file containing all matched functional CDSs.
- Per-sample gene-level RPKM tables.
- A merged gene-level abundance matrix.
- A function-level summed abundance matrix.
- A within-sample relative abundance matrix.
- A column-wise CLR-transformed function matrix.

---

## 2. Required input files

Prepare the following files before running the workflow.

### 2.1 HMM model directory

Download all **27 BileActome `.hmm` files** from the repository directory:

```text
https://github.com/1362996609zby/Bile_acid_v2/tree/main/database%20file
```

Place all downloaded `.hmm` files in one local directory, for example:

```text
hmm_mod/
├── 12ahsdh_final.hmm
├── 12bhsdh_final.hmm
├── 3ahsdh_nadh_final.hmm
├── 3ahsdh_nadph_final.hmm
├── 3ahsdh_steroid_final.hmm
├── 3bhsdh_nadh_final.hmm
├── 3bhsdh_nadph_final.hmm
├── 7ahsdh_nadh_final.hmm
├── 7ahsdh_nadh_nonsteroidal_final.hmm
├── 7ahsdh_nadph_final.hmm
├── 7bhsdh_final.hmm
├── BSH_g_final.hmm
├── BSH_t_final.hmm
├── BaiA1_final.hmm
├── BaiA2_final.hmm
├── BaiB_final.hmm
├── BaiCD_final.hmm
├── BaiE_final.hmm
├── BaiF_final.hmm
├── BaiG_final.hmm
├── BaiH_final.hmm
├── BaiI_final.hmm
├── BaiJ_final.hmm
├── BaiK_final.hmm
├── BaiN_final.hmm
├── BaiO_final.hmm
└── BaiP_final.hmm
```

Each `.hmm` file represents one BileActome functional family. **The HMM filename without the `.hmm` suffix is used as the default function name in all downstream outputs**, including `gene_to_function.tsv`, `abundance_merged_with_function.xlsx`, `FunctionSum_matrix`, `Relative_abundance_matrix`, and `CLR_matrix`.

Example:

```text
BaiCD_final.hmm  →  function = BaiCD_final
BSH_g_final.hmm  →  function = BSH_g_final
7ahsdh_nadh_final.hmm  →  function = 7ahsdh_nadh_final
```

Do not rename the HMM files unless you also want the downstream function names to change.

### 2.2 Non-redundant protein catalog

A FASTA file containing predicted protein sequences from the non-redundant gene catalog.

Example:

```text
Unigenes.protein.cdhit.fa
```

The sequence IDs in this file are used by `hmmscan`.

### 2.3 Non-redundant CDS catalog

A FASTA file containing nucleotide coding sequences corresponding to the protein catalog.

Example:

```text
Unigenes.CDS.cdhit.fa
```

The CDS IDs should match the protein IDs used in the HMM annotation step. If protein IDs and CDS IDs are not identical, use the `--id-regex` option in Step 2 to standardize IDs.

### 2.4 Clean paired-end metagenomic reads

A directory containing paired-end clean reads.

Supported default naming patterns include:

```text
sample.fq1.gz / sample.fq2.gz
sample_R1.clean.fq.gz / sample_R2.clean.fq.gz
sample_R1.fastq.gz / sample_R2.fastq.gz
sample_1.fq.gz / sample_2.fq.gz
```

Users can provide custom R1 glob patterns through Step 3.

---

## 3. Required software and conda environment

### 3.1 Main dependencies

The workflow requires:

- HMMER
- Bowtie2
- Python 3
- Biopython
- pandas
- numpy
- openpyxl

### 3.2 Recommended conda environment

Create a conda environment:

```bash
conda create -n bileactome_env -c conda-forge -c bioconda \
  python=3.10 \
  hmmer \
  bowtie2 \
  biopython \
  pandas \
  numpy \
  openpyxl
```

Activate the environment:

```bash
conda activate bileactome_env
```

Check tools:

```bash
hmmscan -h
hmmpress -h
bowtie2 --version
bowtie2-build --version
python -c "import Bio, pandas, numpy, openpyxl; print('Python dependencies OK')"
```

---

## 4. Directory structure

A recommended project layout is:

```text
project/
├── gene_catalog/
│   ├── Unigenes.protein.cdhit.fa
│   └── Unigenes.CDS.cdhit.fa
├── hmm_mod/
│   ├── BSH_g_final.hmm
│   ├── BSH_t_final.hmm
│   ├── BaiA1_final.hmm
│   ├── BaiCD_final.hmm
│   └── ...  # all 27 downloaded BileActome HMM files
├── clean_reads/
│   ├── sample1.fq1.gz
│   ├── sample1.fq2.gz
│   ├── sample2.fq1.gz
│   └── sample2.fq2.gz
├── hmm_results/
├── bowtie2_map/
├── abundance/
└── results/
```

The scripts can be stored as:

```text
scripts/
├── step1_hmm_annotation.sh
├── step2_extract_matched_cds.py
├── step3_bowtie2_mapping.sh
├── step4_abundance_calc.py
├── step5_merge_abundance_with_function.py
└── step6_function_abundance_transform.py
```

---

## 5. Step-by-step workflow

## Step 0. Prepare HMM indices with hmmpress

### Purpose

`hmmscan` requires pressed HMM databases. Each downloaded BileActome `.hmm` file should be indexed with `hmmpress`, producing `.h3f`, `.h3i`, `.h3m`, and `.h3p` files.

Although Step 1 can automatically run `hmmpress` if index files are missing, users may also prepare them manually.

### Command

```bash
# Put all 27 downloaded BileActome HMM files into one directory first.
for hmm in /path/to/hmm_mod/*.hmm; do
  hmmpress -f "$hmm"
done
```

### Expected outputs

For each HMM file:

```text
BaiCD_final.hmm
BaiCD_final.hmm.h3f
BaiCD_final.hmm.h3i
BaiCD_final.hmm.h3m
BaiCD_final.hmm.h3p
```

The pressed index filenames retain the original HMM basename. Therefore, the family name used downstream remains `BaiCD_final`.

---

## Step 1. HMM annotation of non-redundant proteins

### Script

```text
scripts/step1_hmm_annotation.sh
```

### Purpose

This step runs `hmmscan` to identify proteins in the non-redundant protein catalog that match the custom HMM profiles.

For each HMM family, the script produces a `domtblout` result file and a standard output log.

### Principle

For each `.hmm` file:

```text
protein catalog  →  hmmscan against one HMM  →  matched proteins
```

The output files are named using the HMM filename. For example:

```text
BaiCD_final.hmm  →  BaiCD_final_hmmscan.txt
```

The function name used downstream is therefore `BaiCD_final`.

### Example command

```bash
bash scripts/step1_hmm_annotation.sh \
  --protein-fasta /path/to/gene_catalog/Unigenes.protein.cdhit.fa \
  --hmm-dir /path/to/hmm_mod \
  --out-dir /path/to/hmm_results \
  --cpu 40 \
  --evalue 1e-3 \
  --dom-evalue 1e-3
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--protein-fasta` | Non-redundant protein FASTA file | Required |
| `--hmm-dir` | Directory containing `.hmm` files | Required |
| `--out-dir` | Output directory for hmmscan results | Required |
| `--cpu` | Number of threads | `8` |
| `--evalue` | Sequence-level E-value cutoff | `1e-3` |
| `--dom-evalue` | Domain-level E-value cutoff | `1e-3` |
| `--hmm-pattern` | HMM filename pattern | `*.hmm` |
| `--force` | Force rerun and overwrite previous results | Off |
| `--no-auto-press` | Do not automatically run `hmmpress` | Off |

### Expected outputs

```text
hmm_results/
├── BSH_g_final_hmmscan.txt
├── BSH_g_final_hmmscan.stdout.txt
├── BSH_g_final.done
├── BaiCD_final_hmmscan.txt
├── BaiCD_final_hmmscan.stdout.txt
├── BaiCD_final.done
└── ...
```

### Notes

- The script supports restart/resume.
- If `*.done` and a non-empty result file exist, that HMM family is skipped.
- Temporary files are used to avoid treating incomplete outputs as completed results.

---

## Step 2. Extract matched CDS sequences

### Script

```text
scripts/step2_extract_matched_cds.py
```

### Purpose

This step extracts nucleotide CDS sequences corresponding to proteins that were matched in Step 1.

It also creates a gene-to-function mapping table.

### Principle

The script parses `*_hmmscan.txt` files from Step 1, extracts matched protein IDs, matches these IDs to the CDS FASTA file, and writes:

```text
matched gene IDs + function labels  →  target_cds.fa + gene_to_function.tsv
```

### Example command

```bash
python scripts/step2_extract_matched_cds.py \
  --hmm-result-dir /path/to/hmm_results \
  --cds-fasta /path/to/gene_catalog/Unigenes.CDS.cdhit.fa \
  --out-cds /path/to/target_cds.fa \
  --out-map /path/to/gene_to_function.tsv \
  --force
```

### Example with ID normalization

If protein and CDS IDs require trimming or standardization:

```bash
python scripts/step2_extract_matched_cds.py \
  --hmm-result-dir /path/to/hmm_results \
  --cds-fasta /path/to/gene_catalog/Unigenes.CDS.cdhit.fa \
  --out-cds /path/to/target_cds.fa \
  --out-map /path/to/gene_to_function.tsv \
  --id-regex '^([^|]+)' \
  --force
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--hmm-result-dir` | Directory containing Step 1 hmmscan results | Required |
| `--cds-fasta` | CDS FASTA corresponding to the protein catalog | Required |
| `--out-cds` | Output FASTA of matched functional CDSs | Required |
| `--out-map` | Output gene-to-function TSV | Required |
| `--hmmer-program` | Result format type: `hmmscan` or `hmmsearch` | `hmmscan` |
| `--id-regex` | Regex for standardizing sequence IDs | None |
| `--force` | Overwrite existing outputs | Off |

### Expected outputs

```text
target_cds.fa
gene_to_function.tsv
gene_to_function.tsv.missing_ids.txt
```

`gene_to_function.tsv` format:

```text
gene_id    function
gene0001   BaiCD_final
gene0002   BSH_g_final
gene0003   7ahsdh_nadh_final
```

### Notes

- If many matched IDs are missing from the CDS FASTA, check whether protein IDs and CDS IDs have different formats.
- The `.missing_ids.txt` file helps diagnose ID mismatches.

---

## Step 3. Map clean reads to matched functional CDSs

### Script

```text
scripts/step3_bowtie2_mapping.sh
```

### Purpose

This step maps clean metagenomic reads to `target_cds.fa` using Bowtie2.

### Principle

The extracted target CDSs are used as the reference. Mapping reads back to these CDSs provides read counts for each functional gene.

```text
clean reads  →  Bowtie2 mapping to target_cds.fa  →  per-sample SAM files
```

### Example command

```bash
bash scripts/step3_bowtie2_mapping.sh \
  --target-cds /path/to/target_cds.fa \
  --reads-dir /path/to/clean_reads \
  --out-dir /path/to/bowtie2_map \
  --threads 16
```

### Example with explicit index prefix

```bash
bash scripts/step3_bowtie2_mapping.sh \
  --target-cds /path/to/target_cds.fa \
  --reads-dir /path/to/clean_reads \
  --out-dir /path/to/bowtie2_map \
  --index-prefix /path/to/target_cds_index \
  --threads 16
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--target-cds` | FASTA generated in Step 2 | Required |
| `--reads-dir` | Directory containing paired-end clean reads | Required |
| `--out-dir` | Output directory for mapping results | Required |
| `--index-prefix` | Bowtie2 index prefix | Same directory as `target_cds.fa` |
| `--threads` | Threads for Bowtie2 | `8` |
| `--r1-globs` | R1 file patterns | Common R1 patterns |
| `--force` | Rerun all samples and overwrite outputs | Off |
| `--keep-unal` | Keep unmapped reads in SAM | Off |

### Expected outputs

```text
bowtie2_map/
├── sample1.sam
├── sample1.bowtie2.log
├── sample1.done
├── sample2.sam
├── sample2.bowtie2.log
├── sample2.done
└── ...
```

Bowtie2 index files are also generated:

```text
target_cds_index.1.bt2
target_cds_index.2.bt2
target_cds_index.3.bt2
target_cds_index.4.bt2
target_cds_index.rev.1.bt2
target_cds_index.rev.2.bt2
```

### Notes

- The script supports restart/resume. Existing completed samples are skipped.
- The default behavior uses `--no-unal` to avoid generating very large SAM files.
- If Bowtie2 is killed with signal 9, reduce `--threads`, request more memory, or ensure that `--no-unal` is enabled.

---

## Step 4. Calculate gene-level RPKM abundance

### Script

```text
scripts/step4_abundance_calc.py
```

### Purpose

This step calculates gene-level RPKM abundance for each sample using Bowtie2 SAM files.

### Principle

For each sample:

1. Count mapped reads for each gene.
2. Calculate gene length from `target_cds.fa`.
3. Normalize read counts by gene length and total mapped reads.

The RPKM formula used is:

```text
RPKM = mapped_reads_to_gene / (gene_length_kb × total_mapped_reads_million)
```

Equivalently:

```text
RPKM = read_count × 1e6 / (total_mapped_reads × gene_length_kb)
```

### Example command

```bash
python scripts/step4_abundance_calc.py \
  --sam-dir /path/to/bowtie2_map \
  --cds-fasta /path/to/target_cds.fa \
  --out-dir /path/to/abundance \
  --gene-length-file /path/to/target_cds_length.tsv
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--sam-dir` | Directory containing Step 3 SAM files | Required |
| `--cds-fasta` | `target_cds.fa` from Step 2 | Required |
| `--out-dir` | Output directory for per-sample abundance files | Required |
| `--gene-length-file` | Output gene length table | `target_cds_length.tsv` |
| `--sam-pattern` | SAM filename pattern | `*.sam` |
| `--sample-strip` | String to remove from sample names | `.nohost` |
| `--include-secondary` | Include secondary/supplementary alignments | Off |
| `--id-regex` | Regex for standardizing CDS IDs | None |
| `--force` | Recalculate and overwrite outputs | Off |

### Expected outputs

```text
target_cds_length.tsv
abundance/
├── sample1_rpkm.tsv
├── sample1.done
├── sample2_rpkm.tsv
├── sample2.done
└── ...
```

Each `*_rpkm.tsv` file contains:

```text
Gene    sample1
gene0001    12.45
gene0002    0.00
gene0003    5.78
```

### Notes

- By default, secondary and supplementary alignments are skipped to reduce repeated counting.
- If the same read maps to multiple genes, Bowtie2 alignment behavior and SAM output should be considered when interpreting counts.
- The denominator is total mapped reads to the target CDS reference for that sample.

---

## Step 5. Merge gene-level abundance with function annotation

### Script

```text
scripts/step5_merge_abundance_with_function.py
```

### Purpose

This step merges all per-sample gene-level RPKM files and attaches function labels from `gene_to_function.tsv`.

### Principle

```text
sample-level RPKM files + gene_to_function.tsv
→ gene-level abundance matrix with function annotation
```

### Example command

```bash
python scripts/step5_merge_abundance_with_function.py \
  --abundance-dir /path/to/abundance \
  --function-map /path/to/gene_to_function.tsv \
  --out-xlsx /path/to/abundance_merged_with_function.xlsx \
  --force
```

### Optional TSV output

```bash
python scripts/step5_merge_abundance_with_function.py \
  --abundance-dir /path/to/abundance \
  --function-map /path/to/gene_to_function.tsv \
  --out-xlsx /path/to/abundance_merged_with_function.xlsx \
  --out-tsv /path/to/abundance_merged_with_function.tsv \
  --force
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--abundance-dir` | Directory containing `*_rpkm.tsv` files | Required |
| `--function-map` | `gene_to_function.tsv` from Step 2 | Required |
| `--out-xlsx` | Output Excel file | Required |
| `--out-tsv` | Optional TSV output | None |
| `--join` | Join mode: `inner`, `left`, `outer` | `inner` |
| `--duplicate-policy` | How to handle genes mapped to multiple functions | `first` |
| `--force` | Overwrite existing outputs | Off |

### Expected outputs

```text
abundance_merged_with_function.xlsx
abundance_merged_with_function.tsv   # optional
```

Example matrix:

```text
gene_id    function    sample1    sample2    sample3
gene0001   BaiCD       12.45      8.21       0.00
gene0002   BSH_G       0.00       3.14       1.23
```

### Notes

- The default `inner` join keeps only genes with functional annotation.
- If one gene is assigned to multiple functions, the default policy keeps the first assignment. Use `--duplicate-policy all` if multi-function assignments should be retained.

---

## Step 6. Generate function-level, relative abundance, and CLR matrices

### Script

```text
scripts/step6_function_abundance_transform.py
```

### Purpose

This step aggregates gene-level abundance into function-level matrices and performs compositional transformation.

### Principle

The input table from Step 5 is:

```text
gene_id | function | sample1 | sample2 | ...
```

Step 6 first sums all genes assigned to the same function:

```text
function abundance = sum(RPKM of all genes assigned to that function)
```

Then it calculates within-sample relative abundance:

```text
relative abundance of function i in sample s =
function_i_RPKM_s / sum(all function_RPKM_s in sample s)
```

Finally, it calculates column-wise CLR:

```text
CLR(x_i,s) = ln(x_i,s) - mean_j[ln(x_j,s)]
```

where `i` is a function and `s` is a sample. CLR is performed within each sample column across functions.

### Example command

```bash
python scripts/step6_function_abundance_transform.py \
  --in-table /path/to/abundance_merged_with_function.xlsx \
  --out-func /path/to/FunctionSum_matrix.xlsx \
  --out-rel /path/to/Relative_abundance_matrix.xlsx \
  --out-clr /path/to/CLR_matrix.xlsx \
  --force
```

### Example using TSV input and TSV output

```bash
python scripts/step6_function_abundance_transform.py \
  --in-table /path/to/abundance_merged_with_function.tsv \
  --out-func /path/to/FunctionSum_matrix.tsv \
  --out-rel /path/to/Relative_abundance_matrix.tsv \
  --out-clr /path/to/CLR_matrix.tsv \
  --force
```

### Main parameters

| Parameter | Description | Default |
|---|---|---|
| `--in-table` | Step 5 output table | Required |
| `--out-func` | Function-level summed abundance output | Required |
| `--out-rel` | Relative abundance output | Required |
| `--out-clr` | CLR output | Required |
| `--sheet` | Excel sheet name or index | `0` |
| `--gene-col` | Gene ID column name | `gene_id` |
| `--function-col` | Function column name | `function` |
| `--sample-cols` | Comma-separated sample column names | Auto-detect |
| `--epsilon` | Pseudocount for CLR | `1e-12` |
| `--force` | Overwrite existing outputs | Off |

### Expected outputs

```text
FunctionSum_matrix.xlsx
Relative_abundance_matrix.xlsx
CLR_matrix.xlsx
```

`FunctionSum_matrix.xlsx`:

```text
function    sample1    sample2    sample3
BaiCD_final       120.31     98.11      0.00
BSH_g_final       20.45      14.20      9.77
```

`Relative_abundance_matrix.xlsx`:

```text
function    sample1    sample2    sample3
BaiCD_final       0.42       0.38       0.00
BSH_g_final       0.07       0.05       0.11
```

Each sample column sums to 1.

`CLR_matrix.xlsx`:

```text
function    sample1    sample2    sample3
BaiCD_final       1.23       0.98       -3.45
BSH_g_final       -0.44      -0.61      0.25
```

### Notes

- Use `FunctionSum_matrix` to evaluate total function-level abundance.
- Use `Relative_abundance_matrix` for intuitive composition visualization.
- Use `CLR_matrix` for compositional statistics such as t-tests, PCA, heatmaps, and correlation analyses.
- CLR values indicate whether a function is relatively enriched or depleted compared with the geometric mean of all functions within the same sample.

---

## 6. Complete example workflow

```bash
conda activate bileactome_env

# Step 0. Press HMM files
for hmm in project/hmm_mod/*.hmm; do
  hmmpress -f "$hmm"
done

# Step 1. HMM annotation
bash scripts/step1_hmm_annotation.sh \
  --protein-fasta project/gene_catalog/Unigenes.protein.cdhit.fa \
  --hmm-dir project/hmm_mod \
  --out-dir project/hmm_results \
  --cpu 40 \
  --evalue 1e-3 \
  --dom-evalue 1e-3

# Step 2. Extract matched CDS
python scripts/step2_extract_matched_cds.py \
  --hmm-result-dir project/hmm_results \
  --cds-fasta project/gene_catalog/Unigenes.CDS.cdhit.fa \
  --out-cds project/target_cds.fa \
  --out-map project/gene_to_function.tsv \
  --force

# Step 3. Map reads
bash scripts/step3_bowtie2_mapping.sh \
  --target-cds project/target_cds.fa \
  --reads-dir project/clean_reads \
  --out-dir project/bowtie2_map \
  --threads 16

# Step 4. Calculate RPKM
python scripts/step4_abundance_calc.py \
  --sam-dir project/bowtie2_map \
  --cds-fasta project/target_cds.fa \
  --out-dir project/abundance \
  --gene-length-file project/target_cds_length.tsv

# Step 5. Merge abundance and function annotation
python scripts/step5_merge_abundance_with_function.py \
  --abundance-dir project/abundance \
  --function-map project/gene_to_function.tsv \
  --out-xlsx project/results/abundance_merged_with_function.xlsx \
  --out-tsv project/results/abundance_merged_with_function.tsv \
  --force

# Step 6. Function-level transformation
python scripts/step6_function_abundance_transform.py \
  --in-table project/results/abundance_merged_with_function.xlsx \
  --out-func project/results/FunctionSum_matrix.xlsx \
  --out-rel project/results/Relative_abundance_matrix.xlsx \
  --out-clr project/results/CLR_matrix.xlsx \
  --force
```

---

## 7. Interpretation of final outputs

### 7.1 Gene-level abundance

`abundance_merged_with_function.xlsx` contains the abundance of each matched functional gene across samples.

### 7.2 Function-level summed abundance

`FunctionSum_matrix.xlsx` contains the summed RPKM of all genes assigned to each function.

### 7.3 Relative abundance

`Relative_abundance_matrix.xlsx` contains within-sample proportions across functions.

### 7.4 CLR matrix

`CLR_matrix.xlsx` is recommended for statistical analyses of compositional data.
---
