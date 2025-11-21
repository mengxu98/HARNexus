rm(list = ls())
gc()

source("code/functions/prepare_env.R")

res_dir <- check_dir("../../data/BrainData/processed/GSE186538/")

log_message("Start loading data...")
object <- readRDS(
  file.path(res_dir, "GSE186538.rds")
)

metadata <- object@meta.data
metadata$Cells <- rownames(metadata)
metadata$Dataset <- "GSE186538"
metadata$Technology <- "10X Genomics"
metadata$Sequence <- "snRNA-seq"
metadata$Sample <- metadata$donor_ID
metadata$Sample_ID <- metadata$project_code
annotation <- c(
  # 1. Astrocytes
  "Astro AQP4 CHRDL1" = "Astrocytes",
  "Astro AQP4 GFAP" = "Astrocytes",

  # 2. Oligodendrocytes
  "Oligo OPALIN LAMA2" = "Oligodendrocytes",
  "Oligo OPALIN SLC5A11" = "Oligodendrocytes",
  "Oligo OPALIN LINC01098" = "Oligodendrocytes",
  "Oligo CPXM2 KANK4" = "Oligodendrocytes",

  # 3. Oligodendrocyte progenitor cells (OPC)
  "OPC PDGFRA EGR1" = "Oligodendrocyte progenitor cells",
  "OPC PDGFRA GRIA4" = "Oligodendrocyte progenitor cells",
  "COP GPR17 ADAM33" = "Oligodendrocyte progenitor cells",

  # 4. Microglia
  "Micro C1QB CD83" = "Microglia",
  "Micro C1QB P2RY12" = "Microglia",

  # 5. Myeloid cells
  "Myeloid LSP1 LYZ" = "Myeloid cells",

  # 6. Endothelial cells
  "Endo CLDN5 VWF" = "Endothelial cells",
  "aEndo DKK2 FBLN5" = "Endothelial cells",

  # 7. Vascular mural cells (Pericytes / Smooth muscle / VLMC)
  "PC CLDN5 ABCC9" = "Vascular mural cells",
  "aSMC ACTA2 TAGLN" = "Vascular mural cells",
  "vSMC ABCC9 P2RY14" = "Vascular mural cells",
  "VLMC COL1A1 COL1A2" = "Vascular mural cells",

  # 8. T cells
  "T SKAP1 CD247" = "T cells",

  # 9. Inhibitory neurons
  "InN LAMP5 CHST9" = "Inhibitory neurons",
  "InN LAMP5 KIT" = "Inhibitory neurons",
  "InN LAMP5 NMBR" = "Inhibitory neurons",
  "InN LHX6 AC008415.1" = "Inhibitory neurons",
  "InN MEIS2 SHISAL2B" = "Inhibitory neurons",
  "InN NR2F2 ANO2" = "Inhibitory neurons",
  "InN NR2F2 DDR2" = "Inhibitory neurons",
  "InN NR2F2 MIR4300HG" = "Inhibitory neurons",
  "InN NR2F2 PTPRK" = "Inhibitory neurons",
  "InN NR2F2 SLC17A8" = "Inhibitory neurons",
  "InN PVALB MEPE" = "Inhibitory neurons",
  "InN PVALB PLCL1" = "Inhibitory neurons",
  "InN PVALB PLEKHH2" = "Inhibitory neurons",
  "InN SST ADAMTS12" = "Inhibitory neurons",
  "InN SST EPB41L4A" = "Inhibitory neurons",
  "InN SST NPY" = "Inhibitory neurons",
  "InN SST OTOF" = "Inhibitory neurons",
  "InN VIP ABI3BP" = "Inhibitory neurons",
  "InN VIP CHRNA2" = "Inhibitory neurons",
  "InN VIP NOX4" = "Inhibitory neurons",
  "InN VIP PENK" = "Inhibitory neurons",
  "InN VIP SCML4" = "Inhibitory neurons",
  "InN VIP SCTR" = "Inhibitory neurons",

  # 10. Excitatory neurons
  "CA1 dorsal GRIK1 GRM3" = "Excitatory neurons",
  "CA1 ventral ACVR1C SYT13" = "Excitatory neurons",
  "CA2 CFAP299 HGF" = "Excitatory neurons",
  "CA3 CFAP299 SYN3" = "Excitatory neurons",
  "DG GC PROX1 PDLIM5" = "Excitatory neurons",
  "DG GC PROX1 SGCZ" = "Excitatory neurons",
  "DG MC ARHGAP24 DLC1" = "Excitatory neurons",
  "EC L2 CUX2 CALB1" = "Excitatory neurons",
  "EC L2 CUX2 IL1RAPL2" = "Excitatory neurons",
  "EC L2 CUX2 LAMA3" = "Excitatory neurons",
  "EC L2 CUX2 PDGFD" = "Excitatory neurons",
  "EC L2 RELN BCL11B" = "Excitatory neurons",
  "EC L2 RELN BMPR1B" = "Excitatory neurons",
  "EC L3 PCP4 ADARB2" = "Excitatory neurons",
  "EC L5 BCL11B ADRA1A" = "Excitatory neurons",
  "EC L5 RORB TLL1" = "Excitatory neurons",
  "EC L5 RORB TPBG" = "Excitatory neurons",
  "EC L6 THEMIS CDH13" = "Excitatory neurons",
  "EC L6 THEMIS RGS12" = "Excitatory neurons",
  "EC L6 TLE4 SULF1" = "Excitatory neurons",
  "EC L6b TLE4 CCN2" = "Excitatory neurons",
  "EC L56 TLE4 NXPH2" = "Excitatory neurons",
  "SUB distal FN1 NTNG1" = "Excitatory neurons",
  "SUB proximal ROBO1 COL5A2" = "Excitatory neurons",
  "SUB proximal ROBO1 SEMA3E" = "Excitatory neurons",
  "CR RELN NDNF" = "Inhibitory neurons",
  "Macro F13A1 COLEC12" = "Microglia"
)

metadata$CellType <- annotation[metadata$original_name]
metadata$CellType <- ifelse(
  is.na(metadata$CellType), metadata$original_name, metadata$CellType
)
metadata$Brain_Region <- metadata$region
metadata$Region <- metadata$subregion
metadata$Age <- metadata$donor_age
metadata$Age <- gsub("yr", "", metadata$Age)
metadata$Sex <- metadata$donor_gender
metadata$Sex <- ifelse(metadata$Sex == "F", "Female", "Male")

column_order <- c(
  "Cells", "Dataset", "Technology", "Sequence", "Sample",
  "Sample_ID", "CellType", "Brain_Region", "Region", "Age", "Sex"
)
metadata <- metadata[, column_order]
metadata <- na.omit(metadata)

counts <- GetAssayData(object, layer = "counts")
counts <- counts[, metadata$Cells]
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

log_message("Save data...")
saveRDS(
  object,
  file.path(res_dir, "GSE186538_processed.rds")
)
