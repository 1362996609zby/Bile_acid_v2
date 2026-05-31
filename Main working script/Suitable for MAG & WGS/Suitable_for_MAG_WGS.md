# BileActome MAG/WGS Annotation Workflow

This workflow annotates bile-acid-metabolism genes in microbial genomes, isolate whole-genome sequences (WGS), or metagenome-assembled genomes (MAGs) using BileActome profile HMMs. It performs three main tasks:

1. translate genome nucleotide FASTA files into protein sequences with Prodigal;
2. search predicted proteins against **one BileActome HMM family at a time** using `hmmsearch`;
3. extract matched protein sequences and summarize hit counts per genome and gene family.

This workflow is intended for **genome-level functional annotation**. It differs from the non-redundant gene catalog workflow, which maps metagenomic reads to annotated CDSs and estimates abundance. The final output of this MAG/WGS workflow is a **genome-by-family hit count matrix**, not an abundance matrix.

---

## 1. Required input files

### 1.1 BileActome HMM profiles

Download all **27 BileActome `.hmm` files** from the GitHub repository directory:

```text
https://github.com/1362996609zby/Bile_acid_v2/tree/main/database%20file
```

Place all downloaded `.hmm` files into one local directory, for example:

```text
project_root/BileActome_HMMs/
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

**Important naming rule:** the HMM filename without `.hmm` is used as the default family name in downstream outputs. For example:

```text
BaiE_final.hmm       → family name = BaiE_final
BSH_g_final.hmm      → family name = BSH_g_final
7bhsdh_final.hmm     → family name = 7bhsdh_final
```

This naming rule should be kept consistently through `hmmsearch`, hit extraction, and final matrix generation.

### 1.2 MAG/WGS genome FASTA files

Prepare nucleotide genome FASTA files in `.fna` format:

```text
genomes/
├── genome_001.fna
├── genome_002.fna
└── genome_003.fna
```

Each file should represent one genome, isolate WGS assembly, or MAG. The file stem, such as `genome_001`, will be used as the genome/bin ID in output files.

### 1.3 Current script behavior: one family per run

The current `run_hmmsearch_batch.py` and `extract_hits_from_fasta.py` scripts are designed to process **one HMM family per run**.

Therefore, to annotate all 27 BileActome families, repeat Step 1 and Step 2 separately for each `.hmm` file. For example:

```text
Run 1:  hmm_model = BileActome_HMMs/BaiE_final.hmm
        output_dir = BaiE_final/hmmsearch

Run 2:  hmm_model = BileActome_HMMs/BaiCD_final.hmm
        output_dir = BaiCD_final/hmmsearch

Run 3:  hmm_model = BileActome_HMMs/BSH_g_final.hmm
        output_dir = BSH_g_final/hmmsearch
```

The folder name should match the HMM basename because `count_hmmsearch_tbls.py` treats the first-level folder name as the gene-family name in the final matrix.

Recommended directory layout after running several families:

```text
project_root/
├── BileActome_HMMs/
│   ├── BaiE_final.hmm
│   ├── BaiCD_final.hmm
│   └── ...
├── genomes/
│   ├── genome_001.fna
│   ├── genome_002.fna
│   └── AA/
│       ├── genome_001.faa
│       └── genome_002.faa
├── BaiE_final/
│   └── hmmsearch/
│       ├── genome_001.tbl
│       ├── genome_002.tbl
│       └── hits/
├── BaiCD_final/
│   └── hmmsearch/
│       ├── genome_001.tbl
│       ├── genome_002.tbl
│       └── hits/
└── BSH_g_final/
    └── hmmsearch/
        ├── genome_001.tbl
        ├── genome_002.tbl
        └── hits/
```

---

## 2. Workflow principle

The workflow uses **profile hidden Markov models (profile HMMs)** to detect homologous protein sequences belonging to known bile-acid-metabolism families.

For each genome/MAG and each HMM family:

1. Prodigal predicts proteins from nucleotide contigs.
2. HMMER `hmmsearch` compares predicted proteins against one BileActome HMM profile.
3. HMMER `tblout` files record candidate hits, including target sequence IDs, full-sequence E-values, and HMM scores.
4. Hits passing the E-value threshold, typically `1e-3`, are extracted as FASTA files.
5. Hit counts are summarized across all genomes and all HMM family folders to produce a genome-by-function matrix.

The final count matrix is useful for:

- presence/absence analysis of BileActome genes in MAG/WGS;
- comparison between BileActome and KEGG annotations;
- functional guild classification;
- Bai cluster completeness analysis;
- WGS vs MAG stratified analysis.

---

## 3. Software environment

### 3.1 Create conda environment

```bash
conda create -n bileactome_hmm -c conda-forge -c bioconda \
  python=3.10 hmmer prodigal biopython pandas openpyxl
```

Activate the environment:

```bash
conda activate bileactome_hmm
```

### 3.2 Required packages and tools

| Tool/package | Purpose |
|---|---|
| HMMER (`hmmsearch`) | Search BileActome HMM profiles against predicted proteins |
| Prodigal | Predict protein sequences from MAG/WGS nucleotide FASTA files |
| Python | Run helper scripts |
| Biopython | Read/write FASTA and extract hit sequences |
| pandas | Build count matrices |
| openpyxl | Write Excel output files |

Check installation:

```bash
hmmsearch -h | head
prodigal -h | head
python -c "import Bio, pandas, openpyxl; print('OK')"
```

---

## 4. Step-by-step usage

## Step 0. Download HMM files

### Purpose

Before running the workflow, download all 27 `.hmm` files from the BileActome GitHub `database file` directory and place them into one directory.

### Important note about `hmmpress`

For this MAG/WGS workflow, `hmmsearch` is run against **one single HMM file at a time**. In this use case, `hmmpress` is **not required**. You only need the original `.hmm` files.

Therefore, the first preparation step is simply:

```bash
mkdir -p project_root/BileActome_HMMs
# Download all 27 *.hmm files into project_root/BileActome_HMMs
ls project_root/BileActome_HMMs/*.hmm
```

Expected files include:

```text
BileActome_HMMs/BaiE_final.hmm
BileActome_HMMs/BaiCD_final.hmm
BileActome_HMMs/BSH_g_final.hmm
...
```

---

## Step 1. Predict proteins and run `hmmsearch` for one family

### Script

```bash
python run_hmmsearch_batch.py
```

### Purpose

This step performs two operations:

1. run Prodigal on each `.fna` genome file to generate predicted proteins (`.faa`);
2. run `hmmsearch` between **one specified BileActome HMM profile** and each predicted protein file;
3. produce one `.tbl` result file per genome for that HMM family.

### Principle

```text
genome_001.fna  →  Prodigal  →  genome_001.faa
BaiE_final.hmm + genome_001.faa  →  hmmsearch  →  BaiE_final/hmmsearch/genome_001.tbl
```

### Required edits inside the script

The current script uses path variables inside the Python file. Before each run, edit these variables:

```python
hmm_model = "BileActome_HMMs/BaiE_final.hmm"
genome_dir = "genomes"
protein_dir = os.path.join(genome_dir, "AA")
output_dir = "BaiE_final/hmmsearch"
```

The `output_dir` should use the same basename as the HMM family, without `.hmm`:

```text
BaiE_final.hmm       → output_dir = BaiE_final/hmmsearch
BaiCD_final.hmm      → output_dir = BaiCD_final/hmmsearch
BSH_g_final.hmm      → output_dir = BSH_g_final/hmmsearch
```

### Run example for BaiE_final

```bash
conda activate bileactome_hmm
python run_hmmsearch_batch.py
```

Expected output:

```text
genomes/AA/genome_001.faa
BaiE_final/hmmsearch/genome_001.tbl
BaiE_final/hmmsearch/genome_002.tbl
```

Main commands used internally:

```bash
prodigal -i genome.fna -a genome.faa -p meta -q
hmmsearch --tblout genome.tbl --noali BileActome_HMMs/BaiE_final.hmm genome.faa
```

### Run all 27 families

Because the script processes one HMM family at a time, repeat Step 1 for every `.hmm` file. A simple manual strategy is:

1. edit `hmm_model` and `output_dir` for one family;
2. run `python run_hmmsearch_batch.py`;
3. repeat for the next family.

For example:

```text
BileActome_HMMs/BaiE_final.hmm       → BaiE_final/hmmsearch
BileActome_HMMs/BaiCD_final.hmm      → BaiCD_final/hmmsearch
BileActome_HMMs/BSH_g_final.hmm      → BSH_g_final/hmmsearch
```

### Notes

- `-p meta` is appropriate for MAGs, metagenomic contigs, or unannotated microbial assemblies.
- `--tblout` gives compact tabular output suitable for downstream parsing.
- `--noali` suppresses alignment blocks and keeps output smaller.
- The current `run_hmmsearch_batch.py` script does not require `hmmpress`.

---

## Step 2. Extract matched protein sequences for one family

### Script

```bash
python extract_hits_from_fasta.py
```

### Purpose

This step extracts protein sequences passing the HMMER E-value cutoff for one HMM family.

It performs four operations:

1. read HMMER `.tbl` files from one family-specific `hmmsearch` directory;
2. filter hits by full-sequence E-value;
3. extract matched protein sequences from the corresponding `.faa` files;
4. write hit FASTA files for each genome.

### Required edits inside the script

Before each run, edit these variables so that they match the family just processed in Step 1:

```python
tbl_dir = "BaiE_final/hmmsearch"
genome_dir = "genomes/AA"
output_dir = os.path.join(tbl_dir, "hits")
evalue_threshold = 1e-3
```

For another family:

```python
tbl_dir = "BaiCD_final/hmmsearch"
genome_dir = "genomes/AA"
output_dir = os.path.join(tbl_dir, "hits")
```

### Run example for BaiE_final

```bash
conda activate bileactome_hmm
python extract_hits_from_fasta.py
```

Expected output:

```text
BaiE_final/hmmsearch/hits/genome_001_hits.fasta
BaiE_final/hmmsearch/hits/genome_002_hits.fasta
```

Filtering rule:

```text
full-sequence E-value <= 1e-3
```

### Run all 27 families

Repeat Step 2 for each family-specific `hmmsearch` directory generated in Step 1. The family folder names should remain consistent:

```text
BaiE_final/hmmsearch       → BaiE_final/hmmsearch/hits/
BaiCD_final/hmmsearch      → BaiCD_final/hmmsearch/hits/
BSH_g_final/hmmsearch      → BSH_g_final/hmmsearch/hits/
```

---

## Step 3. Count HMM hits across genomes and gene families

### Script

```bash
python count_hmmsearch_tbls.py
```

### Purpose

This step summarizes all HMMER `.tbl` files across all family folders and generates a genome-by-family count matrix.

It performs four operations:

1. recursively scan all `.tbl` files under the current working directory;
2. treat the first-level folder as the gene family name;
3. count hits passing the E-value threshold for each genome and each gene family;
4. output a genome-by-function count matrix.

### Important parameters inside the script

```python
root_dir = os.getcwd()
EVALUE_THRESHOLD = 1e-3
```

### Required directory structure before running

Run this script from `project_root`, where each family has its own first-level folder:

```text
project_root/
├── BaiE_final/
│   └── hmmsearch/
│       ├── genome_001.tbl
│       └── genome_002.tbl
├── BaiCD_final/
│   └── hmmsearch/
│       ├── genome_001.tbl
│       └── genome_002.tbl
├── BSH_g_final/
│   └── hmmsearch/
│       ├── genome_001.tbl
│       └── genome_002.tbl
└── count_hmmsearch_tbls.py
```

The first-level folder names become the matrix column names:

```text
BaiE_final/	  → column = BaiE_final
BaiCD_final/  → column = BaiCD_final
BSH_g_final/  → column = BSH_g_final
```

### Run

```bash
conda activate bileactome_hmm
cd project_root
python count_hmmsearch_tbls.py
```

Expected output:

```text
hmmsearch_counts_E1e-03.xlsx
```

Output matrix format:

```text
            BaiE_final    BaiCD_final    BSH_g_final    ...
genome_001            2              1              0
genome_002            0              3              1
genome_003            1              0              2
```

Rows are genome/MAG IDs. Columns are BileActome gene families. Values are the number of significant HMM hits.

---

## 5. Final outputs

| Output | Generated by | Description |
|---|---|---|
| `genomes/AA/*.faa` | Step 1 | Predicted proteins from MAG/WGS nucleotide FASTA files |
| `<family>/hmmsearch/*.tbl` | Step 1 | HMMER `tblout` results for each genome and each HMM family |
| `<family>/hmmsearch/hits/*_hits.fasta` | Step 2 | Protein sequences passing the E-value cutoff |
| `hmmsearch_counts_E1e-03.xlsx` | Step 3 | Genome-by-family hit count matrix |

---

## 6. Interpretation of output

### 6.1 Hit FASTA files

The extracted hit FASTA files can be used for:

- manual inspection;
- multiple sequence alignment;
- phylogenetic analysis;
- conserved-domain checking;
- downstream structure modeling or BileActome-X candidate analysis.

### 6.2 Hit count matrix

The count matrix can be converted to:

- presence/absence matrix: `count > 0`;
- Bai cluster completeness score;
- guild classification;
- WGS vs MAG comparison;
- BileActome vs KEGG annotation comparison.

For genome-level pathway completeness, a simple binary transformation is often more robust than raw hit counts:

```text
presence = 1 if hit_count > 0 else 0
```

---

## 7. Recommended downstream analyses

After obtaining the count matrix, common analyses include:

1. **Family-level detection rate**  
   Percentage of genomes carrying each BileActome family.

2. **Bai completeness**  
   Number or fraction of core Bai families detected in each genome.

3. **Guild assignment**  
   Classify genomes into BSH-only, HSDH-only, partial Bai, complete Bai, or mixed guilds.

4. **WGS/MAG stratification**  
   Analyze isolate WGS and MAG representatives separately because MAG contiguity can affect gene-neighborhood interpretation.

5. **Gene-neighborhood analysis**  
   Use genomic coordinates from Prodigal GFF or annotated FASTA headers to test whether Bai/HSDH/BSH genes are colocated, co-oriented, or modularly distributed.

---

## 8. Important notes and limitations

1. **Download all 27 HMM files first.**  
   The workflow assumes that all 27 BileActome `.hmm` files have been downloaded into one local directory.

2. **Do not run `hmmpress` for this workflow unless you only want to validate HMM files.**  
   The current scripts use `hmmsearch` with one HMM file at a time, so pressed HMM index files are not required.

3. **One HMM family is processed at a time.**  
   The first and second scripts need to be run separately for each family by editing `hmm_model`, `output_dir`, and `tbl_dir`.

4. **Family names come from folder names and HMM basenames.**  
   Keep folder names consistent with HMM basenames, for example `BaiE_final.hmm` should use folder `BaiE_final/`.

5. **MAG neighborhood inference requires caution.**  
   MAG fragmentation can split true gene clusters across contigs, causing underestimation of physical proximity.

6. **Hit count is not abundance.**  
   This workflow counts gene hits per genome. It does not estimate metagenomic abundance. For abundance estimation, use the non-redundant gene catalog workflow with Bowtie2 mapping and RPKM/relative abundance/CLR transformation.

7. **E-value cutoff should be reported.**  
   The recommended cutoff here is `1e-3`, but users may adjust it depending on model quality, family specificity, and validation strategy.

---

## 9. Minimal end-to-end example

```bash
# 1. Activate environment
conda activate bileactome_hmm

# 2. Download all 27 HMM files into one directory
mkdir -p BileActome_HMMs
# Manually download all *.hmm files from:
# https://github.com/1362996609zby/Bile_acid_v2/tree/main/database%20file

# 3. Example: run BaiE_final
# Edit run_hmmsearch_batch.py:
# hmm_model = "BileActome_HMMs/BaiE_final.hmm"
# genome_dir = "genomes"
# protein_dir = os.path.join(genome_dir, "AA")
# output_dir = "BaiE_final/hmmsearch"
python run_hmmsearch_batch.py

# 4. Extract BaiE_final hits
# Edit extract_hits_from_fasta.py:
# tbl_dir = "BaiE_final/hmmsearch"
# genome_dir = "genomes/AA"
# output_dir = os.path.join(tbl_dir, "hits")
python extract_hits_from_fasta.py

# 5. Repeat Step 3 and Step 4 for all other 26 HMM families.
#    Keep output folders named after HMM basenames:
#    BaiCD_final/hmmsearch, BSH_g_final/hmmsearch, etc.

# 6. After all families are complete, count all hits
cd project_root
python count_hmmsearch_tbls.py
```

---

