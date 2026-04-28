#### Lake-Zurich-Bacterial-Integrated-Time-Series-LZ-BITS

# Genome architecture predicts proteome allocation in natural bacterial communities

**Cyrill Hofer¹, Lucas Serra Moncadas¹, Jakob Pernthaler¹, Adrian-Stefan Andrei¹***  
¹ Limnological Station, Department of Plant and Microbial Biology, University of Zurich, Kilchberg 8802, Switzerland  
\* Corresponding author: stefan.andrei@limnol.uzh.ch  

[![DOI](https://zenodo.org/badge/1220017213.svg)](https://doi.org/10.5281/zenodo.19735006)

---

## Overview

This repository contains the analysis workflow for the study:

> **Genome architecture predicts proteome allocation in natural bacterial communities**

The project integrates metagenomic and metaproteomic data from a lake time-series to investigate how genome architecture constrains and predicts proteome investment patterns in natural microbial communities.

---

## Repository structure

The repository is organized to reflect the full analysis pipeline:
Data/ → Raw data processing and preparation scripts
Datasets/ → Integrated datasets (combined across data types)
DataAnalysis_SpringBloom2021/ → Analysis scripts for main and supplementary figures


---

## Data workflow

### Data preparation
- Raw data were processed within the `Data/` directory  
- For each data type, curated datasets were generated  
- These were consolidated into unified datasets in the `Datasets/` directory  

### Analysis
- The `DataAnalysis_*` directories contain scripts used to generate:
  - main figures  
  - supplementary figures  
- Each figure has a corresponding analysis workflow  

---

## Data availability

Due to GitHub size limitations, large datasets (e.g. integrated omics tables) are not stored in this repository.

All data are publicly available via Zenodo:

👉 https://doi.org/10.5281/zenodo.19735006

This repository is directly linked to the Zenodo archive via DOI.

---

## Reproducibility

This repository provides:
- analysis scripts  
- data processing workflows  
- figure generation pipelines  

Together with the Zenodo data archive, this enables full reproducibility of the study.

---

## Notes

- Large raw and intermediate data files are intentionally excluded  
- The repository focuses on reproducible workflows and code  
- Data processing is modular and organized by data type  

---

## Contact

**Adrian-Stefan Andrei**  
📧 stefan.andrei@limnol.uzh.ch  

**Cyrill Hofer**
📧 cyrill.hofer@uzh.ch  
---

