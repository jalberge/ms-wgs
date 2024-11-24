setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

sample.sig.genes.table <- read_tsv("../data/MUTSIG2CV/20230509_rsp_mmrf_exomes_jco_su2c_genomes/sample_sig_gene_table.txt")
sig.genes.table <- read_tsv("../data/MUTSIG2CV/20230509_rsp_mmrf_exomes_jco_su2c_genomes/sig_genes.txt")


smm.analysis <- read_tsv(file = '../data/MUTSIG2CV/20230702_noMM_summary_saturation.txt') 
  # filter(n_pts==1034)

final.long.analysis <- read_tsv(file = '../data/MUTSIG2CV/20230701_long_saturation.txt') |> 
  filter(n_pts==1034)


sig.genes.table.main <- sig.genes.table %>%
  filter(npat>10 & q<0.05)

# TODO
# filter mutation artefacts
# filter non expressed genes
# add back forcecalled
# review all mutations

sig.genes.table.main %>% write_tsv("../annot/drivers_20230509_mutsig2cv.txt")

# long.mut.matrix <- sample.sig.genes.table %>%
#   mutate_all(as.character) %>%
#   pivot_longer(cols=-gene, names_to = "mutect_patient_id", values_to = "mutated_or_not")
