# simply merge hg38 SU2C and COMMPASS
library(tidyverse)

#20230509_all_svs_hg38_liftedOver_fixed_ReferencePair.tsv
sv1 <- read_tsv("../data/20230614_MMRF_SV_filtered_results_baseline_hg38_n_776.tsv", num_threads = 4)
sv2 <- read_tsv("../data/20230509_all_svs_hg38_liftedOver_fixed_ReferencePair.tsv", num_threads = 4)

sv.su2c.commpass.final.hg38 <- rbindlist(list(sv1, sv2), fill = TRUE) %>%
  select(-keep)

sv.su2c.commpass.final.hg38.clean <- sv.su2c.commpass.final.hg38 %>% ungroup() %>% filter(stdev1!=0 & stdev2!=0)
sv.su2c.commpass.final.hg38.dirty <- sv.su2c.commpass.final.hg38 %>% filter(!(stdev1!=0 & stdev2!=0))

sv.su2c.commpass.final.hg38 %>% 
  write_tsv("../data/20230509_MMRF_SU2C_MARK_OW_SV_filtered_results_baseline_hg38.tsv")

sv.su2c.commpass.final.hg38.clean %>% 
  write_tsv("../data/20230509_MMRF_SU2C_MARK_OW_SV_filtered_results_baseline_hg38_sd_gt_0.tsv")

sv.su2c.commpass.final.hg38.clean <- read_tsv("../data/20230509_MMRF_SU2C_MARK_OW_SV_filtered_results_baseline_hg38_sd_gt_0.tsv")

# Quick analysis for fun --------------------------------------------------

svs.1 <- sv.su2c.commpass.final.hg38.clean %>% mutate(id=row_number()) %>% transmute(id, individual, chr=chr1, pos=pos1, gene=gene1, class)
svs.2 <- sv.su2c.commpass.final.hg38.clean %>% mutate(id=row_number()) %>% transmute(id, individual, chr=chr2, pos=pos2, gene=gene2, class)

svs.long <- rbindlist(list(svs.1, svs.2))

svs.gene.cluster <- svs.long %>% 
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 

svs.gene.cluster <- svs.long %>% 
  filter(class=="tandem_dup") %>% # DIAPH2
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 

svs.gene.cluster <- svs.long %>% 
  filter(class=="deletion") %>% # PTPRD
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 

svs.gene.cluster <- svs.long %>% 
  # filter(class=="inversion") %>% # EPHA7
  # filter(class=="long_range") %>% # EPHA7 and ZFP36L1
  # filter(class=="inter_chr") %>% # CITED2
  # filter(class=="deletion") %>% # PTPRD
  # filter(class=="tandem_dup") %>% # DIAPH2
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 

svs.gene.cluster <- svs.long %>% 
  # filter(class=="inversion") %>% # EPHA7
  # filter(class=="long_range") %>% # EPHA7 and ZFP36L1
  # filter(class=="inter_chr") %>% # CITED2
  # filter(class=="deletion") %>% # PTPRD
  # filter(class=="tandem_dup") %>% # DIAPH2
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 

svs.gene.cluster <- svs.long %>% 
  # filter(class=="inversion") %>% # EPHA7
  # filter(class=="long_range") %>% # EPHA7 and ZFP36L1
  # filter(class=="inter_chr") %>% # CITED2
  # filter(class=="deletion") %>% # PTPRD
  # filter(class=="tandem_dup") %>% # DIAPH2
  group_by(individual, gene) %>% 
  slice_head(n=1) %>%
  group_by(gene) %>%
  summarise(n=n(), chr=first(chr), min.pos=min(pos), max.pos=max(pos), pts=paste0(individual, collapse = '; ')) %>%
  filter(n>=3) %>%
  arrange(-n) 
