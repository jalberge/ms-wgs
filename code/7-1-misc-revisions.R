library(tidyverse)

source("0_annotate_samples.R")

clinical.and.terra.ref.pairs.samples |>
  filter(Cohort %in% c("CTC", "MARK", "SU2C")) |>
  ungroup() |>
  summarize(median_cov=median(quick_cov_estimate),
            min_cov=min(quick_cov_estimate),
            max_cov=max(quick_cov_estimate))

clinical.and.terra.ref.pairs.samples |>
  filter(Cohort %in% c("SU2C", "MARK")) |>
  ungroup() |>
  summarize(median_cov=median(quick_cov_estimate),
            min_cov=min(quick_cov_estimate),
            max_cov=max(quick_cov_estimate))

clinical.and.terra.ref.pairs.samples |>
  filter(Cohort %in% c("SU2C", "MARK")) |>
  select(Cohort, PANGEA_ID, Tissue, CTF_ID, SCI_ID, Participant_ID, HDP_Participant_ID, HDP_Reference_Pair, Reference_Pair, quick_cov_estimate, quick_cov_estimate_control) |>
  View()
