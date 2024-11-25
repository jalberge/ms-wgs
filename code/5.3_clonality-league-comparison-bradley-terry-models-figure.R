setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# see Andrea's github for code and data
# https://github.com/andrea-poletti-unibo/MS_TiMMing/

library(tidyverse)
library(RColorBrewer)

# data --------------------------------------------------------------------

tm <- read_tsv("../Andrea/data_diff-timing_vs_diff-frequency_SMM_CoMMpass.txt")

head(tm)

tm <- tm |>
  mutate(variable=str_replace(variable, "amp_chr_", "+")) %>%
  mutate(variable=str_replace(variable, "del_chr_", "-")) %>%
  mutate(variable=str_replace(variable, "HyperDiploidy", "HRD")) %>%
  mutate(variable=str_replace(variable, "_MAX", " (MAX)")) %>%
  mutate(variable=str_replace(variable, "_MCL1", " (MCL1)")) %>%
  mutate(variable=str_replace(variable, "_EVI5", " (EVI5)")) %>%
  mutate(variable=str_replace(variable, "_CYLD", " (CYLD)")) %>%
  mutate(variable=str_replace(variable, "_TENT5C", " (TENT5C)")) %>%
  mutate(variable=str_replace(variable, "_XBP1", " (XBP1)")) %>%
  mutate(variable=str_replace(variable, "_MYC_amp", " (MYC)")) %>%
  mutate(variable=str_replace(variable, "_NSD2", " (NSD2)")) %>%
  filter(NDMM>0.05 & SMM > 0.05) |>
  mutate(cn=fct_reorder(variable, mean_timing_com, .desc = TRUE))

g.tm <- ggplot(tm, aes(y=cn, color=)) +
  geom_pointrange(aes(x=mean_timing_com, 
                      xmin=mean_timing_com-sd_com, 
                      xmax=mean_timing_com+sd_com), 
                  color="#333333", fill="#FFFFFF", shape=1) +
  geom_pointrange(aes(x=mean_timing_smm, 
                      xmin=mean_timing_smm-sd_smm, 
                      xmax=mean_timing_smm+sd_smm), 
                  color="#048B9A", alpha=0.7) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank()) +
  labs(x="BT score", y="")

g.tm  

# change HRD for hyperdiploidy
ggsave(filename = "../figures/andrea_timming.png", plot = g.tm, width = 3.2, height = 3.8)
ggsave(filename = "../figures/andrea_timming.pdf", plot = g.tm, width = 3.2, height = 3.8)
