setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(data.table)

# Eventually this is not really used
# We rely on MutSig2CV

# DRIVER SELECTION --------------------------------------------------------

# get drivers of SU2C+CoMMpass+JCO
drivers <- read_tsv("../data/_MUTSIG2CV/20230210.su2c.mark.jco.rsp.mmrf.maf.sig_genes.txt") %>% filter(q<0.1) %>% mutate(dataset="su2c_jco_commpass_exomes")
drivers.commpassonly <- read_tsv("../data/_MUTSIG2CV/MMRF_TvsN_correct_CoveredNormals_NoSerial.aggregated_sig_genes.txt") %>% filter(q<0.1) %>% mutate(dataset="commpass_exomes")
drivers.genomes.only <- read_tsv("../data/_MUTSIG2CV/20230119_su2c_179_commpass_856_hg19_seg_genes.txt") %>% filter(q<0.1) %>% mutate(dataset="su2c_commpass_genomes")
drivers.su2c <- read_tsv("../data/_MUTSIG2CV/20230120_170_patient_su2c_oben_wtc_mark_ctc_no_mmhg19_sig_genes.txt") %>% filter(q<0.1) %>% mutate(dataset="su2c_mgus_smm_only")
drivers.rsp.exomes <- read_tsv("../annot/rsp_commpass_exomes_sig_genes.txt")%>% filter(q<0.1) %>% mutate(dataset="rsp_commpass_exomes")

all.drivers.from.rsp <- read_tsv("../annot/Full_Driver_Table.txt")

drivers.su2c %>% filter(gene %nin% subset(all.drivers.from.rsp, Dataset!="RSP_only")$Gene) %>% View()
drivers %>% filter(gene %nin% all.drivers.from.rsp$Gene) %>% View()
drivers.genomes.only %>% filter(gene %nin% subset(all.drivers.from.rsp, Dataset!="RSP_only")$Gene) %>% View()
drivers.genomes.only %>% filter(gene %nin% all.drivers.from.rsp$Gene) %>% View()

sig.genes <- list(drivers, drivers.commpassonly, drivers.genomes.only, drivers.su2c, drivers.rsp.exomes)

consensus.drivers <- rbindlist(sig.genes, fill=TRUE) %>%
  pivot_wider(id_cols = gene, names_from = dataset, values_from = c(rank, q)) %>%
  left_join(all.drivers.from.rsp, by=c("gene"="Gene"))

clipr::write_clip(consensus.drivers)

consensus.drivers <- consensus.drivers %>% 
  mutate(driver.annotation=case_when(Dataset %in% c("Maura_Walker", "Maura", "Walker") ~ "known driver",
                                     Dataset %in% c("", NA, "RSP_only") & rank_su2c_commpass_genomes > 0 ~ "novel driver genome",
                                     Dataset %in% c("", NA, "RSP_only") & rank_rsp_commpass_exomes > 0 ~ "novel driver exome",
                                     TRUE~"drop"))

exome.union.genome.drivers <- consensus.drivers %>% filter(!is.na(rank_su2c_commpass_genomes)|!is.na(rank_rsp_commpass_exomes))
exome.union.genome.drivers.bar <- consensus.drivers %>% filter(is.na(rank_su2c_commpass_genomes)&is.na(rank_rsp_commpass_exomes))

# 

final.long.analysis <- read_tsv(file = '../data/MUTSIG2CV/20230701_long_saturation.txt')

# manual filter:
# IGLL5 because it's variable region of Ig
filter.manual <- c("IGLL5")

final.long.analysis <- final.long.analysis |> 
  filter(n_pts==1034) |>
  mutate(prevalence=npat/n_pts) |>
  filter(gene %nin% filter.manual) |>
  select(-index)

final.long.analysis |> filter(prevalence>0.01) |> write_tsv("../annot/20230710_mutsig2cv_more_than_1pct.txt")
final.long.analysis |> filter(prevalence<=0.01) |> write_tsv("../annot/20230710_mutsig2cv_less_than_1pct.txt")

# Manual review for synonyms here -----------------------------------------

# genes.of.interest.table %>%
  final.long.analysis %>%
  filter(prevalence<=0.01) %>%
  left_join(all.drivers.from.rsp, by=c("gene"="Gene")) %>%
  relocate(Dataset, Hit) %>% 
  arrange(-prevalence) %>%
  # filter(is.na(Dataset)|Dataset=="RSP_only")
  clipr::write_clip()

# # # manual review secondary names for drivers
# maf <- read_tsv("../maf-for-merge-ctc-su2c/20230203_su2c_oben_wtc_ctc_wm_prm_n_182_patient_level_annotated_hg19.tsv", num_threads = 4, col_select = c(1:22, Protein_Change, t_alt_count, t_ref_count, ccf_hat, ccf_CI_low, ccf_CI_high ))
# rsp.maf <- read_tsv("../data/MMRF_TvsN_correct_CoveredNormals_NoSerial.aggregated.tsv", num_threads = 4)
# # 
# maf |>
#   filter(Variant_Classification %in% non.synonymous) |>
#   filter(Hugo_Symbol %in% c("MALT1")) |>
#   View()
# 
# rsp.maf |>
#   filter(Variant_Classification %in% non.synonymous) |>
#   filter(Hugo_Symbol %in%  c("DNMT3A", "TET2", "ASXL1")) |>
#   select(1:5, Protein_Change, tumor_sample_barcode, ccf_CI95_low, ccf_hat, ccf_CI95_high) |>
#   View()

# # 
# hugo_symbols <- unique(sort(maf$Hugo_Symbol))
# rsp_hugo_symbols <- unique(sort(rsp.maf$Hugo_Symbol))
# 
# exome.union.genome.drivers %>% filter(!(gene %in% hugo_symbols )) %>% pull(gene)
# exome.union.genome.drivers %>% filter(!(gene %in% rsp_hugo_symbols )) %>% pull(gene)

# RNASEH2C => find it in coMMpass
# HIST1H1B => H1-5
# HIST1H1E => H1-4
# RBM16 => SCAF8
# KIAA0182 => GSE1
# MLL4 => KMT2B/2D check from coordinates and number of patients. KMT2D from whole gene name.

synonyms <- exome.union.genome.drivers %>% 
  # filter( gene %in% c("HIST1H1B", "HIST1H1E", "RBM16", "KIAA0182", "MLL4", "WHSC1", "TENT5C") ) %>% 
  mutate(gene=case_when(gene=="HIST1H1B"~"H1-5",
                        gene=="HIST1H1E"~"H1-4",
                        gene=="RBM16"~"SCAF8",
                        gene=="KIAA0182"~"GSE1",
                        gene=="MLL4"~"KMT2D",
                        gene=="WHSC1"~"NSD2", # should not be there but as a security
                        gene=="TENT5C"~"FAM46C", # should not be there but as a security
                        TRUE~gene))

final.drivers.with.synonyms <- synonyms

final.drivers.with.synonyms %>%
  write_tsv("../annot/exome.union.genome.drivers.w.synonyms.tsv")

#"RPL10"    "RBM16"    "DDX3X"    "KIAA0182" "HUWE1"    "IRAK1"    "NKAP"     "MLL4" 
# 
# rsp.maf %>% rowwise %>% 
#   filter(str_detect(paste(HGNC_Previous_Symbols, HGNC_PreviousSymbols1, HGNC_Synonyms, dbNSFP_genename, collapse=" "), "RNASEH2C")) %>% 
#   View
#KAT5/RNASEH2C 	
#11:65479740 (KAT5 IN SU2C)
# => no problemo

