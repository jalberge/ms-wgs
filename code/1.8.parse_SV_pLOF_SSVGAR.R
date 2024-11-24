
library(ggplot2)
library(tidyverse)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# setwd("~/Dropbox (Partners HealthCare)/Projects/ms-wgs/code")

# v7
# sv.m <- read_tsv("../data/JB_MM_SU2C_SV_CoMut_like_plot_input_data_LOF_v7_model_thresh_of_proximity_1000000_min_8_patients.tsv")
# v11
sv.m <- read_tsv("../data/XAVI/JB_MM_SU2C_SV_CoMut_like_plot_input_data_LOF_v11_model_thresh_of_proximity_1000000_min_10_patients.tsv")
sv.m <- sv.m %>% dplyr::rename("Gene"="...1")

patient.order <- colnames(sv.m[-1])
gene.order <- sv.m %>% pull(Gene)

sv.sv <- sv.m %>% 
  pivot_longer(-Gene, names_to = "Sample ID") %>%
  rename(name=Gene) %>%
  relocate(`Sample ID`)
  # mutate(Gene=factor(name, levels=gene.order), 
         # `Sample ID`=factor(`Sample ID`, levels=patient.order))

# signifsv <- read_tsv("../data/Xavi SV LOF results.tsv")
signifsv <- read_tsv("../data/XAVI/new_hit_list_ssvgar_v11.tsv")
signifsv <- signifsv %>% filter(q_value_LOF<0.1)

# Q VALUE FILTER 0.1
sv.sv <- sv.sv %>% filter(name %in% signifsv$gene)

# Utils -------------------------------------------------------------------

# Small function to re-prioritize SVs when manually merging hotspots within 1Mb

prioritize.sv <- function(x, pLOF="At least one SV present; at least one LoF", noLOF="At least one SV present; no LoF") {
  if(length(x)==1) {
    x
  }
  else if (length(x)>1) {
    if (pLOF %in% x) {
      pLOF
    } else if (noLOF %in% x){
      noLOF
    } else
    {
      NA
    }
  }
}


# Check consistency with future files -------------------------------------

HOMEDIR="../data/fig1_tmp/"
METAs=list.files(HOMEDIR, pattern = ".*meta.tsv", full.names = TRUE)
meta.meta <- rbindlist(lapply(METAs, read_tsv))
wide.meta <- meta.meta %>% pivot_wider(id_cols = `Sample ID`)

ref.sample.ids <- sort(unique(wide.meta$`Sample ID`))
svs.sample.ids <- sort(unique(sv.sv$`Sample ID`))

sv.sv.clean <- sv.sv %>%
  mutate(`Sample ID` = case_when(
    str_starts(`Sample ID`, "CTF") ~ str_extract(`Sample ID`, "CTF[0-9]{1,}"),
    str_starts(`Sample ID`, "IID") ~ str_extract(`Sample ID`, "IID_H19606[1-4]"),
    str_starts(`Sample ID`, "MMRF") ~ str_extract(`Sample ID`, "MMRF_[0-9]{4}"),
    str_starts(`Sample ID`, "PD") ~ str_extract(`Sample ID`, "PD[0-9]{5}"),
    `Sample ID`=="pM11077BM_138pos19pos" ~ "pM11077",
    `Sample ID`=="pM11272_BMPCs" ~ "pM11272",
    TRUE ~ `Sample ID`
  ))

svs.clean.sample.ids <- sort(unique(sv.sv.clean$`Sample ID`))

# 1 CTF to pt id (split _ keep 1)
# 2 IID tp pt id (split _ keep 2)
# 3 MMRF to pt id (split _ keep 2)
# Mbp lost on the way
# PDxxx keep str_extract("PD[0-9]{5}")
# pM11077BM_138pos19pos ->
# pM11272_BMPCs ->
setdiff(ref.sample.ids, svs.sample.ids)
setdiff(svs.sample.ids, ref.sample.ids)
setdiff(svs.clean.sample.ids, ref.sample.ids) # great
setdiff(ref.sample.ids, svs.clean.sample.ids)

# CSL151            PIPELINE FAIL
# CSL202            NO SV
# CTF005_CMMCs_57   NO SV passing reviewing
# CTF056_CMMCs_89   NO SV after liftOver (hg38 IGH segment)
# IID_H196063_T0    No SV passing reviewing (part 1)
# IL0303            NO SV
# MBp01t_a, FAIL
# MBp04t_a,FAIL
# MBp09t_a,FAIL
# MBp25t_a,FAIL
# MBp30t_a,         ALL FAIL
# pM10658           NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)
# pM10845           No SV passing reviewing (part 2)
# pM10977           NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)
# pM4175            NO SV (PIPELINE SUCCESS BUT EMPTY RESULT)

# 15 / 180 W/O SV

# RESCUE FROM HG19: CSL162, CSL172, CTF040_CMMCs_264_T1, 
#                   CTF059_BMPCs, pM11272_BMPCs, CTF049_BMPCs (for CSL211)
#                   CTF060_BMPCs, CTF061_BMPCs, CTF062_CTCs_323

# perfect!


# now we need to merge repetitive SVs, and remove IGLL5
sv.sv.clean.grouped <- sv.sv.clean %>%
  mutate(name =
           case_when(
             name %in% c("SP140", "SP140L") ~ "SP140 / SP140L",
             name %in% c("TRAF3", "AMN") ~ "TRAF3 / AMN",
             name %in% c("PRSS2", "MTRNR2L6") ~ "PRSS2 / MTRNR2L6",
             # name == "TRAF2" ~ "TRAF2 SV", # q < 0.25 but not 0.1
             name == "RB1" ~ "RB1 SV",
             TRUE ~ name)) %>%
  filter(!(name %in% "IGLL5")) %>%
  group_by(`Sample ID`, name) %>%
  reframe(value=prioritize.sv(value))
  # group_by(`Sample ID`, name) %>%
  # filter(n()>1L)
  # distinct()

binary.mtx.final <- sv.sv.clean.grouped |>
  pivot_wider(id_cols = name, values_from = value, names_from = `Sample ID`) |>
  tibble::column_to_rownames("name") |>
  mutate(across(everything(), ~as.numeric(ifelse(is.na(.), 0, 1)))) |>
  as.matrix()

# N patients with 1 driver
sum(colSums(binary.mtx.final) >= 1)
# [1] 300
# over 969
300/969

dim(binary.mtx.final)

sv.sv.clean.grouped %>%
  # write_tsv("../data/FINAL_Oct9_SV_Hotspot_Matrix.tsv")
  write_tsv("../data/FINAL_Jun17_SV_Hotspot_Matrix.tsv")

sv.sv.clean <- sv.sv.clean |>
  filter(name %in% signifsv$gene) |>
  mutate(name=factor(name, levels=rev(unique(signifsv$gene))))

sv.patient.order <- sv.sv.clean |>
  pivot_wider(id_cols = `Sample ID`, names_from = name, values_from = value) |>
  # arrange(IGLL5) |>
  arrange(!!!(rev(unique(signifsv$gene)))) |>
  pull(`Sample ID`)
  # mutate(`Sample ID`=fct_reorder(`Sample ID`, ))

sv.sv.clean <- sv.sv.clean |>
  mutate(`Sample ID`=factor(`Sample ID`, levels=sv.patient.order))

p1 <- sv.sv.clean %>%
  ggplot(aes(`Sample ID`, name, fill=value)) +
  geom_tile() +
  scale_fill_manual(values=c("At least one SV present; at least one LoF"="blue", 
                             "At least one SV present; no LoF"="orange",
                             "NA"="white"), na.value = "white") +
  theme_bw() +
  labs(x="Participant", y="", fill="SV predicted function")+
  theme(legend.position = "top", axis.text.x = element_blank(), panel.grid = element_blank())

p1

# signifsv <- read_tsv("../data/Xavi SV LOF results.tsv")
# signifsv <- read_tsv("../data/XAVI/new_hit_list_ssvgar_v11.tsv")

# all(gene.order %in% signifsv$gene)

p2 <- signifsv %>% 
  mutate(Gene=factor(gene, levels=rev(unique(signifsv$gene)))) %>%
  mutate(Q_score=-10*log10(q_value_LOF)) %>%
  ggplot(aes(q_value_LOF, Gene)) +
  geom_col(fill="darkorange") +
  geom_vline(xintercept = c(0.1, 0.25)) +
  scale_x_log10(labels=scales::scientific) +
  labs(x="q score\n(predicted loss-of-function)") +
  theme_bw() 

library(patchwork)  

sv.driver.plot <- p2 + p1 +theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + plot_layout(widths = c(1, 3))
sv.driver.plot

ggsave("../figures/Jun17_reprocess_SV_drivers.png", sv.driver.plot, width = 18/2.54, height = 8/2.54)
ggsave("../figures/Jun17_reprocess_SV_drivers.pdf", sv.driver.plot, width = 18/2.54, height = 8/2.54)
