# assemble SV outputcs from REBC paper method (Chip) from MinimuMM-seq paper, SU2C and external datasets
# hiding some comments which contain RedCapID before github upload
# this will be used by Xavi for his new SV driver method

library(tidyverse)
library(data.table)

setwd("code")

source("0_annotate_samples.R")

# SV with REBC methods don't have the chr prefix in hg38

# SV part 1 (hg38) --------------------------------------------------------
# 2022
sv.part1 <- read_xlsx("../data/SV_Passed_April_7_annotated.xlsx") %>% filter(keep==1)
# we should exclude non-baseline pairs (_b and c)
sv.part1.baseline <- sv.part1 %>% filter(individual %in% clinical.and.terra.ref.pairs.samples$Reference_Pair)

# SV part 2 (hg38) --------------------------------------------------------

sv.part2 <-  read_xlsx("../data/AllWithSVs20221126true_annotated.xlsx") %>% filter(keep==1)
sv.part2.baseline <- sv.part2 %>% filter(individual %in% clinical.and.terra.ref.pairs.samples$Reference_Pair)

# SV part 3 just a few more from sci 21 flowcell --------------------------
# hg38 again
sv.part3 <-  read_xlsx("../data/SV_filtered_results_sci21_2023_annotated.xlsx") %>% filter(keep==1)
sv.part3.baseline <-  sv.part3 %>% filter(individual %in% clinical.and.terra.ref.pairs.samples$Reference_Pair)

# SV from CTC paper hg19 --------------------------------------------------

sv.part4 <- read_xlsx("../../ms-cmmcs/CMMCsPaper2022_ConsensusSV_filtered_results_whitelist_Sept8.xlsx") %>% filter(keep==1) %>% select(-keep)

# remember CTF014 is actually CTF004 (missed longitudinal follow-up when created the CTF ID)
sv.part4[sv.part4[["individual"]]=="CTF014_CMMCs_465", "individual"] <- "CTF004_CMMCs_60"
sv.part4[sv.part4[["individual"]]=="CTF053_CMMCs_124", "individual"] <- "CTF053_BMPCs"

setwd("../maf-for-merge-ctc-su2c/")
sv.part4.hg38 <- liftOverWrapperSV(SV = sv.part4, prefix = "sv.part4", chain = "hg19ToHg38.over.chain", filter_alt_contigs = TRUE, pre_add_chr_prefix = TRUE)
setwd("../code/")
sv.part4.hg38.baseline <- sv.part4.hg38 %>% filter(individual %in% clinical.and.terra.ref.pairs.samples$Reference_Pair)

any( sv.part4.hg38.baseline$individual %in% c(sv.part1.baseline$individual, sv.part2.baseline$individual, sv.part3.baseline$individual))
# FALSE

sv.part4.hg38 %>% write_tsv("../data/CMMCsPaper2022_ConsensusSV_filtered_results_whitelist_Sept8_hg38LiftedOver_use_secondCTF_sample.tsv")
sv.part4.hg38.baseline %>% write_tsv("../data/CMMCsPaper2022_ConsensusSV_filtered_results_whitelist_Sept8_hg38LiftedOver_use_secondCTF_sample_baseline.tsv")

# SV rescue hg19 ----------------------------------------------------------

sv.part5 <- read_xlsx("../data/SV_filtered_hg19_annotated_2023.xlsx") 
sv.part5[sv.part5[["individual"]]=="CTF049_BMPCs", "individual"] <- "CSL211"
sv.part5[sv.part5[["individual"]]=="CTF061_CTCs_106", "individual"] <- "CTF061_BMPCs" # quick hack to rescue 1 tx

setwd("../maf-for-merge-ctc-su2c/")
sv.part5.hg38 <- liftOverWrapperSV(SV = sv.part5, prefix = "sv.part5", chain = "hg19ToHg38.over.chain", filter_alt_contigs = TRUE, pre_add_chr_prefix = TRUE)
setwd("../code/")
sv.part5.hg38.baseline <- sv.part5.hg38 %>% filter(individual %in% clinical.and.terra.ref.pairs.samples$Reference_Pair)

any( sv.part5.hg38.baseline$individual %in% c(sv.part1.baseline$individual, sv.part2.baseline$individual, sv.part3.baseline$individual, sv.part4.hg38.baseline$individual))
# FALSE

sv.part5$individual[sv.part5$individual %nin% clinical.and.terra.ref.pairs.samples$Reference_Pair]

# ALL HG38 SVs

svs <- rbindlist(list(sv.part1.baseline, sv.part2.baseline, sv.part3.baseline, sv.part4.hg38.baseline, sv.part5.hg38.baseline), fill = TRUE)

# missing samples in reference pairs (shouldn't be any)
sort(unique(svs$individual))[ sort(unique(svs$individual)) %nin% clinical.and.terra.ref.pairs.samples$Reference_Pair ]

# reciprocally
clinical.and.terra.ref.pairs.samples$Reference_Pair[ clinical.and.terra.ref.pairs.samples$Reference_Pair %nin% sort(unique(svs$individual)) ]

# CSL202            NO SV
# CTF005_CMMCs_57   NO SV passing reviewing
# CTF056_CMMCs_89   NO SV after liftOver (hg38 IGH segment)
# IID_H196063_T0    No SV passing reviewing (part 1)
# xxxx           NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)
# xxxx           No SV passing reviewing (part 2)
# xxxx           NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)
# xxxx            NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)

# 15 / 180 W/O SV

# RESCUE FROM HG19: xxxx, xxxx, CTF040_CMMCs_264_T1, 
#                   CTF059_BMPCs, xxxx, CTF049_BMPCs (for CSL211)
#                   CTF060_BMPCs, CTF061_BMPCs, CTF062_CTCs_323

svs %>% write_tsv("../data/20230509_all_svs_hg38_liftedOver_fixed_ReferencePair.tsv")

