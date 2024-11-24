setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)


# Before MutSig2CV Rerun --------------------------------------------------


`%nin%` <- Negate("%in%")

maf.mutsig <- read_tsv("../data/20230503.su2c.mark.jco.rsp.mmrf.hg19.nowm.nomgip.maf", num_threads = 8)

length(unique(maf.mutsig$Tumor_Sample_Barcode))
# 1037
# down to #1030
participant.list <- read_tsv("../data/export_MM-like_score.tsv")

dim(participant.list)

colnames(maf.mutsig)

source("0_annotate_samples.R")
clinical.and.terra.ref.pairs.samples

maf.mutsig <- maf.mutsig |>
  mutate(candidate = case_when( Tumor_Sample_Barcode %in% participant.list$ID ~ Tumor_Sample_Barcode,
                                startsWith(prefix = "MMRF", Tumor_Sample_Barcode) ~ str_extract(Tumor_Sample_Barcode, "MMRF_[0-9]{4}"),
                                startsWith(prefix = "SMM", Tumor_Sample_Barcode) ~ str_replace_all(Tumor_Sample_Barcode, "-", "_"),
                                TRUE ~ Tumor_Sample_Barcode
                                ) )

# which ones are in mutsig2cv but not in the matrix?
maf.mutsig.tumors <- unique(sort(maf.mutsig$candidate))

maf.mutsig.tumors[maf.mutsig.tumors %nin% participant.list$ID]

# 1037
# pM11003 --> removed because AL --> 1036
# pM4175 --> removed because CLL --> 1035
# "MMRF_1889" "MMRF_2363" --> removed because follow up! --> 1033
# "SMM_076_Tumor" "SMM_081_Tumor" "SMM_093_Tumor" --> probably not clean. have to review. --> 1030
# CSL193 in MAF, is SMM_003_Tumor in List (SMM_003_tUmor not in maf)
missing.still <- maf.mutsig.tumors[maf.mutsig.tumors %nin% participant.list$ID]

# conversely 
equi <- data.frame(ID=participant.list$ID[participant.list$ID %nin% maf.mutsig.tumors])
equi.match <- equi |> 
  left_join(clinical.and.terra.ref.pairs.samples |> select(Participant_ID, HDP_Participant_ID), by=c("ID"="Participant_ID"))
equi.match$HDP_Participant_ID %in% maf.mutsig.tumors
# only SMM_TUMOR_003 is in participant list but not in maf

FINALEXTRACLEANMAF <- maf.mutsig |> filter(Tumor_Sample_Barcode %nin% c("pM11003", "pM4175", "MMRF_1889_2_BM", "MMRF_2363_2_BM", "SMM-076-Tumor", "SMM-081-Tumor", "SMM-093-Tumor"))
length(unique(FINALEXTRACLEANMAF$Tumor_Sample_Barcode))
FINALEXTRACLEANMAF$candidate <- NULL

FINALEXTRACLEANMAF |> write_tsv("../data/20240618_export_1030_baseline.maf")
  

# After MutSig2CV rerun ---------------------------------------------------

results <- read_tsv("/Users/jean-baptiste/Dropbox (Partners HealthCare)/2_Projects/ms-wgs/data/_MUTSIG2CV/20240619_1030Run/outdir/sig_genes.txt")

# load old final
final.long.analysis <- read_tsv(file = '../data/_MUTSIG2CV/20230701_long_saturation.txt')
filter.manual <- c("IGLL5")
final.long.analysis <- final.long.analysis |> 
  filter(n_pts==1034) |>
  mutate(prevalence=npat/n_pts) |>
  filter(gene %nin% filter.manual) |>
  select(-index)
# final.long.analysis |> filter(prevalence>0.01) |> write_tsv("../annot/20230710_mutsig2cv_more_than_1pct.txt")
# final.long.analysis |> filter(prevalence<=0.01) |> write_tsv("../annot/20230710_mutsig2cv_less_than_1pct.txt")

original.results <- final.long.analysis
new.results <- results |> filter(q<0.1)

new.results$prevalence = new.results$npat/1030
new.results |> filter( gene %nin% original.results$gene )
#  "EHD1" -> New (not found by Maura or Walker)
# "RFTN1"  -> also found by Maura Walker
# "ZNF292"  -> also found by Maura Walker
# (these are the new hits from last MutSig2CV 1030 participants).
original.results |> filter( gene %nin% new.results$gene )
# IDH1 missing -> make sure it's removed from figure
# BTG1 missing ->  make sure it's removed from figure
# KRT2 missing -> not in figure, make sure it's not in the table
# POT1 missing ->  make sure it's removed from figure


new.results |> filter(prevalence>0.01) |> write_tsv("../annot/20240619_mutsig2cv_more_than_1pct.txt")
new.results |> filter(prevalence<=0.01) |> write_tsv("../annot/20240619_mutsig2cv_less_than_1pct.txt")
