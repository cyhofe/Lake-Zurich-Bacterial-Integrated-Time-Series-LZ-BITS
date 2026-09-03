#### Lake-Zurich-Bacterial-Integrated-Time-Series-LZ-BITS

# Genome architecture predicts proteome allocation in natural bacterial communities

**Cyrill Hofer¹, Lucas Serra Moncadas¹, Jakob Pernthaler¹, Adrian-Stefan Andrei¹**\*  
¹ Limnological Station, Department of Plant and Microbial Biology, University of Zurich, Kilchberg 8802, Switzerland  
\* Corresponding author: stefan.andrei@limnol.uzh.ch

**Archived analysis release and supporting data:**  
👉 DOI:/10.5281/zenodo.19735006

---

## Overview

This repository contains the analysis workflow for the study:

> **Genome architecture predicts proteome allocation in natural bacterial communities**

The project integrates metagenomic, metaproteomic, flow-cytometric, and environmental data from a Lake Zurich time series to investigate how genome architecture constrains and predicts proteome investment patterns in natural microbial communities.

---

## Repository structure

The repository is organized to reflect the full analysis pipeline:

- `Data/` → Data processing, preparation, quality-control, and data-type-specific analysis workflows
- `Datasets/` → Curated and integrated datasets combining information across data types
- `DataAnalysis_SpringBloom2021/` → Analysis scripts for main and supplementary figures

Depending on the analysis, scripts use processed or curated files directly from `Data/`, integrated datasets from `Datasets/`, or a combination of both.

---

## Data workflow

### Data preparation

- Raw data were processed using workflows contained in the `Data/` directory
- Data-type-specific processed and curated outputs were generated within the corresponding `Data/` subdirectories
- Where required, information from multiple data types was consolidated into integrated datasets in the `Datasets/` directory

### Analysis

- The `DataAnalysis_SpringBloom2021/` directory contains scripts used to generate:
  - main figures
  - supplementary figures
  - statistical analyses
- Depending on the specific analysis, input data are taken directly from `Data/`, from integrated datasets in `Datasets/`, or from both
- Each figure and analysis has a corresponding analysis workflow

---

## Data availability

Due to GitHub size limitations, large raw, intermediate, processed, and integrated data files are not stored directly in this repository.

Raw Lake Zurich omics data are deposited in the corresponding domain-specific repositories:

- **Metagenomic sequencing:** ENA accession `PRJEB95916`
- **Metaproteomics:** PRIDE accession `PXD070123`

The archived analysis release and additional supporting datasets are available via Zenodo:

👉 DOI: /10.5281/zenodo.19735006

The Zenodo archive contains the data required by the analysis workflows in this repository, including:

- processed and curated files from `Data/`
- integrated datasets from `Datasets/`
- archived analysis scripts and workflows
- supporting files required for figure generation and statistical analyses
- Lake Zurich flow-cytometry data
- MockD6300 metagenomic and metaproteomic benchmarking data

Raw Lake Zurich flow-cytometry data are provided through Zenodo as `.fcs` files together with the corresponding CytExpert `.xit` files. Measurements were acquired using a **CytoFLEX S flow cytometer (Beckman Coulter)** and analyzed using **CytExpert v2.4.0.28**.

The flow-cytometry archive includes both the high-resolution integrated SpringBloom 2021 measurements corresponding to the 37 LZ-BITS samples and the complementary biweekly depth-resolved Lake Zurich sampling series.

A **ZymoBIOMICS Microbial Community Standard D6300 (MockD6300)** was used as an independent benchmarking dataset for the metagenomic assembly and metaproteomic workflows.

To avoid mixing mock-community benchmarking samples with the environmental Lake Zurich samples in the primary ENA and PRIDE repositories, the MockD6300 data are provided separately through Zenodo.

The MockD6300 data available on Zenodo include:

- paired-end Illumina raw sequencing reads
- Nanopore raw sequencing reads
- raw and corrected draft genome assemblies
- metagenomic assembly benchmarking outputs
- raw proteomics `.mzML` files
- proteomics processing and benchmarking outputs

---

## Data repository

**Lake Zurich sequencing data:**  
ENA – European Nucleotide Archive  
`PRJEB95916`

**Lake Zurich proteomics data:**  
PRIDE – Proteomics Identifications Database  
`PXD070123`

**Archived analysis data, Lake Zurich flow cytometry, MockD6300 benchmarking data, and supporting files:**  
Zenodo  
DOI: /10.5281/zenodo.19735006

---

## Reproducibility

This repository provides:

- analysis scripts
- data processing workflows
- statistical analysis workflows
- figure generation pipelines

Depending on the analysis, these workflows use processed or curated files directly from `Data/`, integrated datasets from `Datasets/`, or both.

The corresponding data are available from the appropriate public repositories:

- ENA for raw Lake Zurich metagenomic sequencing data
- PRIDE for raw Lake Zurich metaproteomic data
- Zenodo for processed and integrated analysis data, Lake Zurich flow-cytometry data, MockD6300 benchmarking data, and supporting files

Together, the GitHub repository and the associated Zenodo, ENA, and PRIDE deposits provide the code and data required to reproduce the analyses presented in the study.

---

## Notes

- Large raw, intermediate, processed, and integrated data files are intentionally excluded from GitHub
- The repository focuses on reproducible workflows and analysis code
- Analyses may use data directly from `Data/`, integrated datasets from `Datasets/`, or both
- Raw environmental Lake Zurich sequencing data are deposited in ENA
- Raw environmental Lake Zurich metaproteomic data are deposited in PRIDE
- Lake Zurich flow-cytometry data are deposited in Zenodo
- MockD6300 sequencing, genome-assembly, and proteomics benchmarking data are deposited in Zenodo
- Data processing is modular and organized by data type

---

## Contact

**Adrian-Stefan Andrei**  
📧 stefan.andrei@limnol.uzh.ch

**Cyrill Hofer**  
📧 cyrill.hofer@uzh.ch
