setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(survival)
library(viridis)
library(survminer)

source("0_annotate_samples.R")

# load score from SMM
smm.scores <- read_tsv("../data/20241110_MM-like_IMWG_cytogenet_scores.tsv")

# survival SU2C
survival.SU2C <- read_tsv("../data/export_MM-like_score_SMM_survival_SU2C.tsv")
survival.SU2C$ID=survival.SU2C$`Sample ID`
# survival MARK
survival.mark <- read_csv("../JB/Modified dataframes/Mark_Signature_SU2C_only.csv")

mark.su2c.survival <- full_join(survival.SU2C, survival.mark)
mark.su2c.survival$study <- ifelse(mark.su2c.survival$cohort%in% c("SU2C", "CTC"), "SU2C", "JCO")

HR.Gain1q <- mark.su2c.survival |>
  pivot_longer(cols = c( starts_with("amp_") , starts_with("t1"), "KRAS", "NRAS", "FAM46C", starts_with("del_")),
               names_to = "mutation", values_to = "status") |>
  nest_by(mutation) |>
  mutate(model=list(coxph(Surv(Days, Progression_status) ~ status, data=data))) |>
  summarise(tidy(model, conf.int=TRUE, exponentiate=TRUE)) |>
  ungroup() |>
  mutate(q.value=p.adjust(p.value, method="fdr")) |>
  relocate(p.value, .before = q.value) |>
  relocate(statistic, .before = p.value)

HR.Gain1q |> write_tsv("../data/Export_Gain1Q_and_others_MM-like-significance.tsv")
