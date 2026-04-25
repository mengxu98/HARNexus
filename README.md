# HARNexus

## Overview

![HARNexus overview](code/overview.svg)

HARNexus is an integrative framework for mapping the regulatory **nexus** between human accelerated regions (HARs) and downstream target genes at single-cell resolution.

HARs are rapidly evolved noncoding regulatory elements linked to human-specific brain traits. HARNexus combines HAR annotations, TF motif information, and sc/snRNA-seq data to build cell type–specific regulatory networks (CSNs) and prioritize HAR target genes.

Using independent Hi-C evidence, HARNexus supports both target-gene prediction accuracy and cell-type specificity. The framework was applied to a large human brain atlas (~2.24M sc/snRNA-seq profiles; 44 regions, 15 developmental stages, 9 cell types) and extended to cross-species analyses with snATAC-seq to study regulatory rewiring in human brain evolution.

Overall, HARNexus provides a practical and scalable way to study HAR function in brain development and evolution.

## Repository

This repository contains analysis pipelines and plotting code.

### Data dependency

The large processed sc/snRNA-seq and snATAC-seq objects are managed in the companion repository:
https://github.com/mengxu98/BrainOmicsData

The scripts expect that repository's processed outputs under `../../data/BrainOmicsData/` relative to this repository. Local input files for HAR annotations, motifs, Hi-C support, and selected gene sets are stored under `data/`.

### Basic usage

Run the numbered scripts from the repository root. The HAR-TF step should be run before network analyses because it creates `results/har_tf/tfs.csv`, which is reused by the cerebellum, atlas, and species workflows.

1. `01_har_tf_pairs_analysis.sh`
2. `02_datasets_processing.sh`
3. `03_datasets_integration.sh`
4. `04_cerebellum_networks.sh`
5. `05_hic.sh`
6. `06_har_csn_atlas.sh`
7. `06_species_networks.sh`
8. `08_plotting.sh`
9. `09_supplementary_tables.sh`
