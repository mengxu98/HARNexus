#!/bin/bash

set -e

# Source the enhanced download utilities and log functions
source "code/functions/utils.sh"

code_dir="code/datasets"
overwrite="${1:-F}"

should_process() {
  local target_file="$1"
  if [[ "$overwrite" =~ ^([Tt]|[Tt][Rr][Uu][Ee]|1)$ ]]; then
    return 0
  fi
  if [ ! -f "$target_file" ]; then
    return 0
  fi
  return 1
}


# GSE103723
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103723
# code:
if should_process "../../data/BrainData/processed/GSE103723/GSE103723_processed.rds"; then
  log_message "Processing GSE103723 data..."
  Rscript $code_dir/GSE103723.R
  log_success "GSE103723 data processed successfully!"
else
  log_message "GSE103723 data already processed!"
fi


# GSE104276
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276
# code:
if should_process "../../data/BrainData/processed/GSE104276/GSE104276_processed.rds"; then
  log_message "Processing GSE104276 data..."
  Rscript $code_dir/GSE104276.R
  log_success "GSE104276 data processed successfully!"
else
  log_message "GSE104276 data already processed!"
fi


# # GSE126836 (no age information)
# # paper: https://doi.org/10.1016/j.cell.2019.05.006
# # pmid: https://pubmed.ncbi.nlm.nih.gov/31178122/
# # data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836
# # code:
# if should_process "../../data/BrainData/processed/GSE126836/GSE126836_processed.rds"; then
#   log_message "Processing GSE126836 data..."
#   Rscript $code_dir/GSE126836.R
#   log_success "GSE126836 data processed successfully!"
# else
#   log_message "GSE126836 data already processed!"
# fi


# GSE186538
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538
# code:
if should_process "../../data/BrainData/processed/GSE186538/GSE186538_processed.rds"; then
  log_message "Processing GSE186538 data..."
  Rscript $code_dir/GSE186538.R
  log_success "GSE186538 data processed successfully!"
else
  log_message "GSE186538 data already processed!"
fi


# GSE199762
# paper: https://doi.org/10.1038/s41586-023-06981-x
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/38122823
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199762
# dbGaP, accession: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003509.v1.p1
# code: https://github.com/massisnascimento/ECstream
if should_process "../../data/BrainData/processed/GSE199762/GSE199762_processed.rds"; then
  log_message "Processing GSE199762 data..."
  Rscript $code_dir/GSE199762.R
  log_success "GSE199762 data processed successfully!"
else
  log_message "GSE199762 data already processed!"
fi


# GSE204683 (Multiome: snRNA-seq + snATAC-seq (GSE204684))
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204683
# ATAC-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204684
# CELLxGENE (RRID: SCR_021059) data (h5ad):
# https://cellxgene.cziscience.com/collections/ceb895f4-ff9f-403a-b7c3-187a9657ac2c
# code: https://doi.org/10.5281/zenodo.7703253
bash $code_dir/download_GSE204683.sh

if should_process "../../data/BrainData/processed/GSE204683/GSE204683_processed.rds"; then
  log_message "Processing GSE204683 data..."
  Rscript $code_dir/GSE204683.R
  log_success "GSE204683 data processed successfully!"
else
  log_message "GSE204683 data already processed!"
fi


# GSE212606
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/37824614
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212606
# code:
if should_process "../../data/BrainData/processed/GSE212606/GSE212606_processed.rds"; then
  log_message "Processing GSE212606 data..."
  Rscript $code_dir/GSE212606.R
  log_success "GSE212606 data processed successfully!"
else
  log_message "GSE212606 data already processed!"
fi


# GSE217511
# paper: https://doi.org/10.1038/s41467-022-34975-2
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36509746
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217511
if should_process "../../data/BrainData/processed/GSE217511/GSE217511_processed.rds"; then
  log_message "Processing GSE217511 data..."
  Rscript $code_dir/GSE217511.R
  log_success "GSE217511 data processed successfully!"
else
  log_message "GSE217511 data already processed!"
fi


# GSE67835
# paper: https://doi.org/10.1073/pnas.1507125112
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/26060301
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835
if should_process "../../data/BrainData/processed/GSE67835/GSE67835_processed.rds"; then
  log_message "Processing GSE67835 data..."
  Rscript $code_dir/GSE67835.R
  log_success "GSE67835 data processed successfully!"
else
  log_message "GSE67835 data already processed!"
fi


# GSE81475
# paper: https://doi.org/10.1016/j.celrep.2016.08.038
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/27568284
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81475
# code:
if should_process "../../data/BrainData/processed/GSE81475/GSE81475_processed.rds"; then
  log_message "Processing GSE81475 data..."
  Rscript $code_dir/GSE81475.R
  log_success "GSE81475 data processed successfully!"
else
  log_message "GSE81475 data already processed!"
fi


# GSE97942 (GSE97887 + GSE97930)
# paper: https://doi.org/10.1126/sciadv.adg3754
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/29227469
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97942
# GSE97887 (scTHS-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97887
# GSE97930 (snDrop-seq): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930
# code:
if should_process "../../data/BrainData/processed/GSE97942/GSE97942_processed.rds"; then
  log_message "Processing GSE97942 data..."
  Rscript $code_dir/GSE97942.R
  log_success "GSE97942 data processed successfully!"
else
  log_message "GSE97942 data already processed!"
fi


# Li_et_al_2018
# paper: https://doi.org/10.1126/science.aat7615
# pmid: https://pubmed.ncbi.nlm.nih.gov/30545854/
# data:
# code:
if should_process "../../data/BrainData/processed/Li_et_al_2018/Li_et_al_2018_processed.rds"; then
  log_message "Processing Li_et_al_2018 data..."
  Rscript $code_dir/Li_et_al_2018.R
  log_success "Li_et_al_2018 data processed successfully!"
else
  log_message "Li_et_al_2018 data already processed!"
fi


# Nowakowski_et_al_2017
# paper: https://doi.org/10.1126/science.aap8809
# pmid: https://pubmed.ncbi.nlm.nih.gov/29217575/
# data:
# code:
if should_process "../../data/BrainData/processed/Nowakowski_et_al_2017/Nowakowski_et_al_2017_processed.rds"; then
  log_message "Processing Nowakowski_et_al_2017 data..."
  Rscript $code_dir/Nowakowski_et_al_2017.R
  log_success "Nowakowski_et_al_2017 data processed successfully!"
else
  log_message "Nowakowski_et_al_2017 data already processed!"
fi


# BTSatlas
# Contains multiple datasets:
# "AllenM1", "EGAD00001006049", "EGAS00001006537",
# "GSE178175", "GSE168408",
# "GSE144136", "GSE202210"
# first run BTSatlas-1.py to create compatible h5ad file
# then run BTSatlas-2.R to split the data into multiple datasets
if [ ! -f "../../data/BrainData/processed/BTSatlas/BTS_atlas_compatible.h5ad" ]; then
  log_message "Processing BTSatlas data..."
  python $code_dir/BTSatlas-1.py
  log_success "BTSatlas data processed successfully!"
else
  log_message "BTSatlas compatible h5ad file already exists!"
fi

log_message "Processing BTSatlas sub-datasets..."
Rscript $code_dir/BTSatlas_2.R
log_success "BTSatlas sub-datasets processed successfully!"

# AllenM1
# paper:
# pmid:
# data:
# code:
if should_process "../../data/BrainData/processed/AllenM1/AllenM1_processed.rds"; then
  log_message "Processing AllenM1 data..."
  Rscript $code_dir/AllenM1.R
  log_success "AllenM1 data processed successfully!"
else
  log_message "AllenM1 data already processed!"
fi


# EGAD00001006049
if should_process "../../data/BrainData/processed/EGAD00001006049/EGAD00001006049_processed.rds"; then
  log_message "Processing EGAD00001006049 data..."
  Rscript $code_dir/EGAD00001006049.R
  log_success "EGAD00001006049 data processed successfully!"
else
  log_message "EGAD00001006049 data already processed!"
fi


# EGAS00001006537
if should_process "../../data/BrainData/processed/EGAS00001006537/EGAS00001006537_processed.rds"; then
  log_message "Processing EGAS00001006537 data..."
  Rscript $code_dir/EGAS00001006537.R
  log_success "EGAS00001006537 data processed successfully!"
else
  log_message "EGAS00001006537 data already processed!"
fi


# GSE178175
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175
# code:
if should_process "../../data/BrainData/processed/GSE178175/GSE178175_processed.rds"; then
  log_message "Processing GSE178175 data..."
  Rscript $code_dir/GSE178175.R
  log_success "GSE178175 data processed successfully!"
else
  log_message "GSE178175 data already processed!"
fi


# GSE168408
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
# code:
if should_process "../../data/BrainData/processed/GSE168408/GSE168408_processed.rds"; then
  log_message "Processing GSE168408 data..."
  Rscript $code_dir/GSE168408.R
  log_success "GSE168408 data processed successfully!"
else
  log_message "GSE168408 data already processed!"
fi


# GSE144136
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
# code:
if should_process "../../data/BrainData/processed/GSE144136/GSE144136_processed.rds"; then
  log_message "Processing GSE144136 data..."
  Rscript $code_dir/GSE144136.R
  log_success "GSE144136 data processed successfully!"
else
  log_message "GSE144136 data already processed!"
fi


# GSE202210
# paper:
# pmid:
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202210
# code:
if should_process "../../data/BrainData/processed/GSE202210/GSE202210_processed.rds"; then
  log_message "Processing GSE202210 data..."
  Rscript $code_dir/GSE202210.R
  log_success "GSE202210 data processed successfully!"
else
  log_message "GSE202210 data already processed!"
fi


# SCR_016152
# SCR_016152 contains GSE207334 and Ma_et_al_2022,
# which are from the same study: https://www.ncbi.nlm.nih.gov/pubmed/36007006
# paper: https://doi.org/10.1126/science.abo7257
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9614553/
# pmid: https://www.ncbi.nlm.nih.gov/pubmed/36007006
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207334
# http://resources.sestanlab.org/PFC/
# https://brainscope.gersteinlab.org/
# code:
# all samples not include GSE207334
# object_all <- readRDS(
#   "../../data/BrainData/raw/GSE207334/PFC_snRNAseq_liftover.rds"
# )

# GSE207334 (Multiome: snRNA-seq + snATAC-seq)
if should_process "../../data/BrainData/processed/GSE207334/GSE207334_processed.rds"; then
  log_message "Processing GSE207334 data..."
  Rscript $code_dir/GSE207334.R
  log_success "GSE207334 data processed successfully!"
else
  log_message "GSE207334 data already processed!"
fi

# Ma_et_al_2022
if should_process "../../data/BrainData/processed/Ma_et_al_2022/Ma_et_al_2022_processed.rds"; then
  log_message "Processing Ma_et_al_2022 data..."
  Rscript $code_dir/Ma_et_al_2022.R
  log_success "Ma_et_al_2022 data processed successfully!"
else
  log_message "Ma_et_al_2022 data already processed!"
fi


# GSE235493 (Multiome: snRNA-seq + snATAC-seq, Macaque)
# paper: https://doi.org/10.1016/j.neuron.2025.04.025
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235493
# code: https://doi.org/10.5281/zenodo.15243470
# note: not found human data in GSE235493

# HYPOMAP
# paper: https://doi.org/10.1038/s41586-024-08504-8
# code:
#   https://github.com/lsteuernagel/HYPOMAP
#   https://github.com/georgiedowsett/HYPOMAP
#   https://github.com/lsteuernagel/scIntegration
#   https://github.com/mrcepid-rap
# data: https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480
# data(spatial): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278848

# SomaMut
# papaer: https://doi.org/10.1038/s41586-025-09435-8
# code: ~
# data: https://publications.wenglab.org/SomaMut/
if should_process "../../data/BrainData/processed/SomaMut/SomaMut_processed.rds"; then
  log_message "Processing SomaMut data..."
  Rscript $code_dir/SomaMut.R
  log_success "SomaMut data processed successfully!"
else
  log_message "SomaMut data already processed!"
fi

# PRJCA015229 (Multiome: snRNA-seq + snATAC-seq, Human + Macaque)
# paper: https://doi.org/10.1016/j.xgen.2024.100703
# code: https://github.com/KIZ-SubLab/ACC-sn-Multiomes
# data: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015229
if should_process "../../data/BrainData/processed/PRJCA015229/PRJCA015229_processed.rds"; then
  log_message "Processing PRJCA015229 data..."
  Rscript $code_dir/PRJCA015229.R
  log_success "PRJCA015229 data processed successfully!"
else
  log_message "PRJCA015229 data already processed!"
fi

# ROSMAP (Religious Order Study (ROS) or the Rush Memory and Aging Project (MAP))
# paper: https://doi.org/10.1016/j.cell.2023.08.039
# code: https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC/
# data: https://compbio.mit.edu/ad_aging_brain/

bash $code_dir/download_rosmap_processed_data.sh
bash $code_dir/download_rosmap_ucsc_snRNAseq.sh
bash $code_dir/download_rosmap_ucsc_snRNAseqsnATACseq.sh
# bash $code_dir/download_rosmap_ucsc_snATACseq.sh
# bash $code_dir/download_rosmap_ucsc_snATACseq_Epigenomic.sh
if should_process "../../data/BrainData/processed/ROSMAP/ROSMAP_processed.rds"; then
  log_message "Processing ROSMAP data..."
  Rscript $code_dir/ROSMAP.R
  log_success "ROSMAP data processed successfully!"
else
  log_message "ROSMAP data already processed!"
fi


# datasets integration
if should_process "../../data/BrainData/processed/integration/integration_list.rds"; then
  log_message "Processing datasets integration..."
  Rscript $code_dir/datasets_integration.R
  log_success "Datasets integration processed successfully!"
else
  log_message "Datasets integration already processed!"
fi
