# BileActome

BileActome is a curated HMM-based resource for annotating microbial bile acid metabolism genes. This repository provides the BileActome HMM database files and example workflow scripts for two common use scenarios: genome-resolved annotation of MAG/WGS datasets and community-level profiling based on non-redundant metagenomic gene catalogs.

## Repository structure

```text
Bile_acid_v2/
├── database file/
├── Main working script/
│   ├── Suitable for MAG & WGS/
│   └── Suitable_for_non_redundant_gene_catalogs/
└── Database building scripts/
```

## Recommended starting point

For practical use of BileActome, please go directly to the **Main working script** folder.

This folder contains two workflow pipelines:

### 1. Suitable for MAG & WGS

Use this workflow if your input data are microbial genomes, including:

- metagenome-assembled genomes (MAGs)
- isolate whole-genome sequences (WGS)
- predicted protein files from genome-resolved datasets

This workflow is designed for genome-resolved annotation of bile acid metabolism genes using the BileActome HMM database.

The required input files, detailed running steps, parameters, and output files are fully described in the README file inside the **Suitable for MAG & WGS** folder.

### 2. Suitable_for_non_redundant_gene_catalogs

Use this workflow if your input data are non-redundant metagenomic gene catalogs.

This workflow is designed for community-level profiling of bile acid metabolism genes based on non-redundant gene sets and metagenomic read mapping.

The required input files, detailed running steps, parameters, and output files are fully described in the README file inside the **Suitable_for_non_redundant_gene_catalogs** folder.

## BileActome database files

Before running either workflow, please download all `.hmm` files in the **database file** folder.

These `.hmm` files are the BileActome database files and are required for annotation.

The detailed usage of these database files is described in the README files inside each workflow folder under **Main working script**.

## Purpose of the workflow scripts

The scripts provided in the **Main working script** folder are intended to:

- provide example parameters for using BileActome;
- demonstrate complete workflows for the two major use cases;
- help users understand how BileActome can be applied to MAG/WGS datasets and non-redundant gene catalogs.

These scripts are not intended to be the only possible way to use BileActome. Users may modify the scripts, adjust parameters, or integrate the BileActome HMM database into their own bioinformatic pipelines according to their research needs.

## Database building scripts

The **Database building scripts** folder contains scripts used by the authors during the construction of the BileActome database.

These scripts are provided for transparency and reproducibility of database construction.

For routine use of BileActome, users do not need to run the scripts in the **Database building scripts** folder.

## General workflow

A typical use of BileActome follows these steps:

1. Download the BileActome `.hmm` files from the **database file** folder.
2. Choose the appropriate workflow under **Main working script**:
   - use **Suitable for MAG & WGS** for genome-resolved MAG/WGS annotation;
   - use **Suitable_for_non_redundant_gene_catalogs** for non-redundant gene catalog profiling.
3. Read the README file inside the selected workflow folder.
4. Prepare the required input files according to the workflow-specific instructions.
5. Run the scripts using the recommended parameters or adapt them to your own dataset.
6. Interpret the output files according to the workflow-specific documentation.

## Notes

BileActome is designed as a reusable annotation resource for microbial bile acid metabolism. The example workflows in this repository are provided to facilitate reproducible use of the database in two common microbiome analysis scenarios, but users are encouraged to adapt the database and scripts to their own data structures and computational environments.
