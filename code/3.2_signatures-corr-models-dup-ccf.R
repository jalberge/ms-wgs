setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggforestplot)
library(tidyverse)
library(ComplexHeatmap)
library(ggpmisc)

# this file contains all correlative analyses from HDP's output
# eventually, needs to be splitted into several smaller files
# could be: 
# 1- dedup mutations / samples
# 2- review mutant vaf (supplementary notes) for artifacts versus novel unmatched signatures
# 2- correlate loadings with disease stage
# 3- extract CCF of signatures
# 4- extract signature profile of early mutations (duplicated)

source("0_annotate_samples.R")

HOMEDIR="../data/fig1_tmp/"
METAs=list.files(HOMEDIR, pattern = ".*meta.tsv", full.names = TRUE)
meta.meta <- rbindlist(lapply(METAs, read_tsv))
meta <- meta.meta %>% pivot_wider(id_cols = `Sample ID`)

meta.sigs <- meta %>% filter(Cohort %nin% "Bustoros")

# these one needs to remain controlled-access
sigs <- read.table("../data/_MM_Sigs_HDP/sig_exposures_final.txt", row.names = 1)
maf.sbs <- read_tsv("../data/_MM_Sigs_HDP/20240304_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv.gz", num_threads = 7)

# Separate again clustered, IG, total, etc. -------------------------------
# HDP's output needs to reprocessing (in particular sample names)

# column names are weird (1 added at the end?)
first.clustered <- which(endsWith(colnames(sigs), "_clustered"))[1]
first.ig <- which(endsWith(colnames(sigs), "_IG"))[1]
first.1 <- which(endsWith(colnames(sigs), "CSL1411"))[1]

sigs.default <- sigs[, 1:(first.clustered-1)]
sigs.clustered <-  sigs[, first.clustered:(first.ig-1)]
sigs.IG <- sigs[, first.ig:(first.1-1)]
sigs.1 <- sigs[, first.1:ncol(sigs)]

# REQUIRES MAF
stats <- tibble(maf.sbs) |> 
  group_by(Participant_ID, SampleID, Tumor_Sample_Barcode, Tumor_Sample_Barcode_combined) |> 
  count() |> 
  ungroup() |>
  mutate(Group = case_when(
        str_detect(SampleID, "_clustered$") ~ "clustered",
        str_detect(SampleID, "_IG$") ~ "IG",
        TRUE ~ "Genome"))

tumor.sample.barcode.combined.2.patient <- stats |> 
  select(Participant_ID, Tumor_Sample_Barcode_combined) |> 
  distinct() |> 
  ungroup()


# Color palette -----------------------------------------------------------

#SBS84-like and N17 are the same
sbs.palette = 
  c("N9"="#333333", "N12"="#5c5c5c", "N16"="#858585", "N18"="#adadad", "N19"="#d6d6d6", "N20"="#dddddd",
  "N13"="#C7E9C0", "N15"="#A1D99B", "SBS8"="#74C476", "SBS16"="#31A354", "SBS18"="#006D2C",
  "SBS1"="#FA8072", "SBS5"="#B22222", "SBS40"="#B22222", "SBS2"="#FFA500", "SBS13"="#EEEE00", 
  "SBS17a"="#EEAEEE", "SBS17b"="#B452CD", "SBS9"="#BFEFFF", "N17"="#00BFFF", "SBS84-like"="#00BFFF", "SBS85"="#1874CD")
sbs.order <- c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15", "SBS8", 
               "SBS16", "SBS18", "SBS1", "SBS5", "SBS40", "SBS2", "SBS13", "SBS17a", 
               "SBS17b", "SBS9", "N17", "SBS84-like", "SBS85")
# combined
combined.samples.and.samples <- sigs.default |>
  rownames_to_column("SBS") |>
  pivot_longer(-SBS, names_to = "sample_id", values_to = "weight") |>
  mutate(SBS=factor(SBS, levels=sbs.order))

ggplot(combined.samples.and.samples, aes(sample_id, weight, fill=SBS)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values=sbs.palette)


# average at patient level ------------------------------------------------

patient.level.sigs <- combined.samples.and.samples |>
  left_join(tumor.sample.barcode.combined.2.patient, by=c("sample_id"="Tumor_Sample_Barcode_combined")) |>
  group_by(Participant_ID, SBS) |>
  summarise(weight=mean(weight))

patient.order.apobec <- patient.level.sigs |> 
  filter(SBS %in% c("N9")) |>
  group_by(Participant_ID) |>
  summarise(weight=sum(weight)) |>
  arrange(-weight)

ggplot(patient.level.sigs, aes(factor(Participant_ID, levels=patient.order.apobec$Participant_ID), weight, fill=SBS)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values=sbs.palette)

# Annotate patients -------------------------------------------------------

meta <- read_tsv("../data/FINAL_meta_columns_matrix.tsv")

meta$`Sample ID`

meta.hdp <- meta %>% 
  filter(Cohort %nin% c("Bustoros")) %>%
  left_join(clinical.and.terra.ref.pairs.samples %>% select(HDP_Participant_ID, HDP_Reference_Pair, Participant_ID), by=c("Sample ID" = "Participant_ID")) %>%
  mutate(HDP_Participant_ID=ifelse(is.na(HDP_Participant_ID), `Sample ID`, HDP_Participant_ID))

sigs.annot <- patient.level.sigs |>
  inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID"))

sigs.annot <- sigs.annot |> 
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA"
  ))

sigs.annot$Stage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"))
sigs.annot$SimpleStage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"), labels=c("MGUS/SMM", "MGUS/SMM", "MM"))

sigs.annot |> write_tsv("../data/_MM_Sigs_HDP/sigs.annot.patients.tsv")

# test without sequencing artefacts ---------------------------------------------

artefacts <- c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15") 

sigs.annot.no.artefact <- sigs.annot |>
  ungroup() |>
  filter(SBS %nin% artefacts) |>
  group_by(Participant_ID) |>
  mutate(weight=weight/sum(weight))


# Group specific signatures -----------------------------------------------
# Gender - Age don't matter (and Female bias with MAF)

# model SBS 
# weight ~ SimpleStage + IMWG + WGS

## Model ----------

stats.sigs.annot <- sigs.annot |> 
  ungroup() |>
  mutate(WGS=factor(Cohort=="MMRF", levels=c(FALSE, TRUE))) |>
  mutate(SimpleStage=factor(SimpleStage, levels=c("MM", "MGUS/SMM"))) |>
  mutate(IMWG=factor(IMWG, levels=c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) |>
  group_by(SBS) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG + WGS , data = .), conf.int=TRUE)) |>
  filter(term!="(Intercept)") |>
  group_by(term) |>
  mutate(p.adj=p.adjust(p.value, method="fdr")) |>
  ungroup()# adjust for each signature tested

stats.sigs.annot <- stats.sigs.annot |>
  mutate(term=case_when(
    term=="SimpleStageMGUS/SMM"~"MGUS/SMM [ref: MM]",
    term=="WGSTRUE"~"MMRF [ref: Standard WGS]",
    startsWith(term, "IMWG") ~ paste0(str_remove(term, "IMWG"), " [ref: Unclassified]"),
    TRUE ~ term
  ))

stats.sigs.annot.noMMRF <- sigs.annot |> 
  filter(Cohort!="MMRF") |>
  group_by(SBS) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG, data = .), conf.int=TRUE))

## Forest plots per disease stage and group --------------------------------

# Draw a forestplot of cross-sectional, linear associations
p.sbs.apobec <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS2", "SBS13")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "APOBEC"
) 

p.sbs.clock <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS1", "SBS5")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Clock-like"
) 

p.sbs.aid.shm <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS9", "SBS85", "SBS84-like", "N17")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "AID / SHM"
)

p.unknown.ros <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS8", "SBS16", "SBS18")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Unknown / ROS"
) 

p.unknown.17 <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS17a", "SBS17b")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "SBS17a/b"
) 

p.new.artefacts <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Novel / Possible artefacts"
)

library(patchwork)

page.sbs.graphs <- p.sbs.apobec +
  p.sbs.clock + 
  p.sbs.aid.shm + 
  p.unknown.ros +
  p.unknown.17 +
  p.new.artefacts +
  plot_layout(ncol = 2) &
  labs(x="Beta [95% CI]", color="") &
  scale_color_manual(values=sbs.palette) &
  scale_x_continuous(limits=c(-0.12, 0.19)) &
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        legend.position = "bottom"
        )

page.sbs.graphs

ggsave("../figures/sbs-models-per-group.pdf", 
       page.sbs.graphs, 
       width = 7, 
       height = 9)

ggsave("../figures/sbs-models-per-group.png", 
       page.sbs.graphs, 
       width = 7, 
       height = 9)

# Genome-wide mutations signatures --------------------------------------------

## Model ----------
stats.sigs.annot <- sigs.annot |> 
  ungroup() |>
  filter(Cohort!="MMRF") |>
  mutate(SimpleStage=factor(SimpleStage, levels=c("MM", "MGUS/SMM"))) |>
  mutate(IMWG=factor(IMWG, levels=c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) |>
  group_by(SBS) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG, data = .), conf.int=TRUE)) |>
  filter(term!="(Intercept)") |>
  group_by(term) |>
  mutate(p.adj=p.adjust(p.value, method="fdr")) |>
  ungroup()# adjust for each signature tested

stats.sigs.annot <- stats.sigs.annot |>
  mutate(term=case_when(
    term=="SimpleStageMGUS/SMM"~"MGUS/SMM [ref: MM]",
    startsWith(term, "IMWG") ~ paste0(str_remove(term, "IMWG"), " [ref: Unclassified]"),
    TRUE ~ term
  ))


## Forest plots per disease stage and group --------------------------------
  

# Draw a forestplot of cross-sectional, linear associations
p.sbs.apobec <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS2", "SBS13")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "APOBEC"
) 

p.sbs.clock <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS1", "SBS5")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Clock-like"
) 

p.sbs.aid.shm <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS9", "SBS85", "SBS84-like", "N17")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "AID / SHM"
)

p.unknown.ros <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS8", "SBS16", "SBS18")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Unknown / ROS"
) 

p.unknown.17 <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS17a", "SBS17b")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "SBS17a/b"
) 

p.new.artefacts <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Novel / Possible artefacts"
)

library(patchwork)

page.sbs.graphs <- p.sbs.apobec +
  p.sbs.clock + 
  p.sbs.aid.shm + 
  p.unknown.ros +
  p.unknown.17 +
  p.new.artefacts +
  plot_layout(ncol = 2) &
  labs(x="Beta [95% CI]", color="") &
  scale_color_manual(values=sbs.palette) &
  scale_x_continuous(limits=c(-0.12, 0.19)) &
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        legend.position = "bottom"
  )

ggsave("../figures/sbs-models-per-group-noMMRF.pdf", 
       page.sbs.graphs, 
       width = 7, 
       height = 9)

ggsave("../figures/sbs-models-per-group-noMMRF.png", 
       page.sbs.graphs, 
       width = 7, 
       height = 9)

# Clustered mutations -----------------------------------------------------

combined.samples.and.samples <- sigs.clustered |>
  rename_all(~str_remove(., "_clustered")) |>
  rownames_to_column("SBS") |>
  pivot_longer(-SBS, names_to = "sample_id", values_to = "weight") |>
  mutate(SBS=factor(SBS, levels=sbs.order))

patient.level.sigs <- combined.samples.and.samples |>
  left_join(tumor.sample.barcode.combined.2.patient, by=c("sample_id"="Tumor_Sample_Barcode_combined")) |>
  group_by(Participant_ID, SBS) |>
  summarise(weight=mean(weight))

sigs.annot <- patient.level.sigs |>
  inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID"))

sigs.annot <- sigs.annot |> 
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA"
  ))

sigs.annot$Stage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"))
sigs.annot$SimpleStage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"), labels=c("MGUS/SMM", "MGUS/SMM", "MM"))

sigs.annot |> write_tsv("../data/_MM_Sigs_HDP/sigs.annot.patients.clustered.tsv")


## Model ----------

stats.sigs.annot <- sigs.annot |> 
  ungroup() |>
  mutate(WGS=factor(Cohort=="MMRF", levels=c(FALSE, TRUE))) |>
  mutate(SimpleStage=factor(SimpleStage, levels=c("MM", "MGUS/SMM"))) |>
  mutate(IMWG=factor(IMWG, levels=c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) |>
  group_by(SBS) |>
  # do(tidy(lm(weight ~ SimpleStage + IMWG + WGS + Gender + Age, data = .), conf.int=TRUE)) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG + WGS , data = .), conf.int=TRUE)) |>
  filter(term!="(Intercept)") |>
  group_by(term) |>
  mutate(p.adj=p.adjust(p.value, method="fdr")) |>
  ungroup()# adjust for each signature tested

stats.sigs.annot <- stats.sigs.annot |>
  mutate(term=case_when(
    term=="SimpleStageMGUS/SMM"~"MGUS/SMM [ref: MM]",
    term=="WGSTRUE"~"MMRF [ref: Standard WGS]",
    startsWith(term, "IMWG") ~ paste0(str_remove(term, "IMWG"), " [ref: Unclassified]"),
    TRUE ~ term
  ))

stats.sigs.annot.noMMRF <- sigs.annot |> 
  filter(Cohort!="MMRF") |>
  group_by(SBS) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG, data = .), conf.int=TRUE))


## Forest plots per disease stage and group --------------------------------

# Draw a forestplot of cross-sectional, linear associations
p.sbs.apobec <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS2", "SBS13")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "APOBEC"
) 

p.sbs.clock <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS1", "SBS5")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Clock-like"
) 

p.sbs.aid.shm <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS9", "SBS85", "SBS84-like", "N17")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "AID / SHM"
)

p.unknown.ros <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS8", "SBS16", "SBS18")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Unknown / ROS"
) 

p.unknown.17 <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS17a", "SBS17b")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "SBS17a/b"
) 

p.new.artefacts <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Novel / Possible artefacts"
)

library(patchwork)

page.sbs.graphs.clustered <- p.sbs.apobec +
  p.sbs.clock + 
  p.sbs.aid.shm + 
  p.unknown.ros +
  p.unknown.17 +
  p.new.artefacts +
  plot_layout(ncol = 2) &
  labs(x="Beta [95% CI]", color="") &
  scale_color_manual(values=sbs.palette) &
  scale_x_continuous(limits=c(-0.16, 0.21)) &
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        legend.position = "bottom"
  )

page.sbs.graphs.clustered

ggsave("../figures/sbs-models-per-group.clustered.pdf", 
       page.sbs.graphs.clustered, 
       width = 7, 
       height = 9)

ggsave("../figures/sbs-models-per-group.clustered.png", 
       page.sbs.graphs.clustered, 
       width = 7, 
       height = 9)



# IG mutations -----------------------------------------------------

combined.samples.and.samples <- sigs.IG |>
  rename_all(~str_remove(., "_IG")) |>
  rownames_to_column("SBS") |>
  pivot_longer(-SBS, names_to = "sample_id", values_to = "weight") |>
  mutate(SBS=factor(SBS, levels=sbs.order))

patient.level.sigs <- combined.samples.and.samples |>
  left_join(tumor.sample.barcode.combined.2.patient, by=c("sample_id"="Tumor_Sample_Barcode_combined")) |>
  group_by(Participant_ID, SBS) |>
  summarise(weight=mean(weight))

sigs.annot <- patient.level.sigs |>
  inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID"))

sigs.annot <- sigs.annot |> 
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA"
  ))

sigs.annot$Stage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"))
sigs.annot$SimpleStage <- factor(sigs.annot$Stage, levels=c("MGUS", "SMM", "MM"), labels=c("MGUS/SMM", "MGUS/SMM", "MM"))

sigs.annot |> write_tsv("../data/_MM_Sigs_HDP/sigs.annot.patients.IG.tsv")


## Model ----------

stats.sigs.annot <- sigs.annot |> 
  ungroup() |>
  mutate(WGS=factor(Cohort=="MMRF", levels=c(FALSE, TRUE))) |>
  mutate(SimpleStage=factor(SimpleStage, levels=c("MM", "MGUS/SMM"))) |>
  mutate(IMWG=factor(IMWG, levels=c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) |>
  group_by(SBS) |>
  # do(tidy(lm(weight ~ SimpleStage + IMWG + WGS + Gender + Age, data = .), conf.int=TRUE)) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG + WGS , data = .), conf.int=TRUE)) |>
  filter(term!="(Intercept)") |>
  group_by(term) |>
  mutate(p.adj=p.adjust(p.value, method="fdr")) |>
  ungroup()# adjust for each signature tested

stats.sigs.annot <- stats.sigs.annot |>
  mutate(term=case_when(
    term=="SimpleStageMGUS/SMM"~"MGUS/SMM [ref: MM]",
    term=="WGSTRUE"~"MMRF [ref: Standard WGS]",
    startsWith(term, "IMWG") ~ paste0(str_remove(term, "IMWG"), " [ref: Unclassified]"),
    TRUE ~ term
  ))

stats.sigs.annot.noMMRF <- sigs.annot |> 
  filter(Cohort!="MMRF") |>
  group_by(SBS) |>
  do(tidy(lm(weight ~ SimpleStage + IMWG, data = .), conf.int=TRUE))


## Forest plots per disease stage and group --------------------------------

# Draw a forestplot of cross-sectional, linear associations
p.sbs.apobec <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS2", "SBS13")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "APOBEC"
) 

p.sbs.clock <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS1", "SBS5")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Clock-like"
) 

p.sbs.aid.shm <- forestplot(
  df = stats.sigs.annot |> filter(SBS %in% c("SBS9", "SBS85", "SBS84-like", "N17")),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "AID / SHM"
)

p.unknown.ros <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS8", "SBS16", "SBS18")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Unknown / ROS"
) 

p.unknown.17 <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("SBS17a", "SBS17b")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "SBS17a/b"
) 

p.new.artefacts <- forestplot(
  df = stats.sigs.annot |> 
    filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15")),
  # group_by(term) |>
  # mutate(label=paste0(term, "\n", expression(beta), "=", format(estimate, 2))) |>
  # ungroup(),
  name = term,
  estimate = estimate,
  pvalue = p.adj,
  psignif = 0.1,
  se = std.error,
  colour = SBS,
  title = "Novel / Possible artefacts"
)

library(patchwork)

page.sbs.graphs.ig <- p.sbs.apobec +
  p.sbs.clock + 
  p.sbs.aid.shm + 
  p.unknown.ros +
  p.unknown.17 +
  p.new.artefacts +
  plot_layout(ncol = 2) &
  labs(x="Beta [95% CI]", color="") &
  scale_color_manual(values=sbs.palette) &
  scale_x_continuous(limits=c(-0.35, 0.145)) &
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        legend.position = "bottom"
  )

page.sbs.graphs.ig

ggsave("../figures/sbs-models-per-group.ig.pdf", 
       page.sbs.graphs.ig, 
       width = 7, 
       height = 9)

ggsave("../figures/sbs-models-per-group.ig.png", 
       page.sbs.graphs.ig, 
       width = 7, 
       height = 9)

# Re-do SBS plot from Tim ---------------------------------------------------

## dedup and average -------------------------------------------------------

dedup.maf.sbs <- tibble(maf.sbs) |>
  # filter(Participant_ID=="CTF001") |>
  group_by(Participant_ID, Chromosome, Start_position) |>
  slice_head(n=1) |>
  ungroup() |>
  mutate(Group = case_when(
    str_detect(SampleID, "_clustered$") ~ "clustered",
    str_detect(SampleID, "_IG$") ~ "IG",
    TRUE ~ "Genome")) |>
  inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID")) |>
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) |>
  mutate(SimpleStage=factor(Stage, levels=c("MGUS", "SMM", "MM"), labels=c("MGUS/SMM", "MGUS/SMM", "MM"))) |>
  mutate(Guess=case_when(
    SBS9 + SBS85 + N17 > 0.85  ~ "AID/SHM",
    SBS2 + SBS13 > 0.85 ~ "APOBEC",
    SBS1 + SBS5 + SBS40 > 0.85 ~ "Clock-like",
    TRUE ~ "All<0.85")) |>
  mutate(clonal_landau = ccf_median>=0.85 & !is.na(ccf_median), non_clonal_landau= ccf_median<0.85 & !is.na(ccf_median))

# threshold inspection
ggplot( dedup.maf.sbs, aes(SBS2)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS13)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS1)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS5)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS40)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS9)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(SBS85)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))
ggplot( dedup.maf.sbs, aes(N17)) + geom_histogram() + scale_y_log10() + geom_vline(xintercept = c(0.5, 0.85, 0.9))


dedup.maf.sbs |> write_tsv("../data/_MM_Sigs_HDP/20240415_deduplicated_maf.tsv")

# dedup.maf.sbs <- read_tsv("../data/_MM_Sigs_HDP/20240415_deduplicated_maf.tsv")

sbs.cols <- c("N9", "N12", 
          "N13", "N15", "N16", "N17", "N18", "N19", "N20", "SBS1", "SBS13", 
          "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS2", "SBS40", "SBS5", 
          "SBS8", "SBS85", "SBS9")


# VAF/maf for unmatched signatures ------------------------------------------
# 2024 06 13

dedup.maf.sbs.with.probable.process <- 
  dedup.maf.sbs[ apply(dedup.maf.sbs[, sbs.cols], 1, function(x) { any(x>0.5)}), ]

max.sbs <- apply(dedup.maf.sbs.with.probable.process[, sbs.cols], 1, which.max)
dedup.maf.sbs.with.probable.process$which_max <- max.sbs
dedup.maf.sbs.with.probable.process$which_max_SBS = sbs.cols[max.sbs]

dedup.maf.sbs.with.probable.process.ccf_hat <- dedup.maf.sbs.with.probable.process |>
  # sample_frac(size=0.1) |>
  mutate(which_max_SBS=fct_recode(which_max_SBS, `SBS84-like`="N17")) |>
  mutate(which_max_SBS_group=case_when(
    which_max_SBS %in% c("N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    which_max_SBS %in% c("N9", "N13", "N15") ~ "Unmatched",  
    which_max_SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    which_max_SBS %in% c("SBS18") ~ "ROS",
    which_max_SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    which_max_SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    which_max_SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    which_max_SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA"
  )) |>
  mutate(maxi_group=case_when(which_max_SBS_group %in% c("AID/SHM", "APOBEC", "Clock-like", "Unmatched") ~ "1",
                              TRUE ~ "2"))

ccf.1.line <- ggplot(dedup.maf.sbs.with.probable.process.ccf_hat |> filter(maxi_group=="1"),
                     aes(which_max_SBS, ccf_hat, fill=which_max_SBS)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale="width") +
  scale_fill_manual(values=sbs.palette) +
  facet_grid(~ which_max_SBS_group, scales = "free_x", space = "free_x") +
  labs(x="", y="Cancer Cell Fraction Maximum Estimate") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ccf.1.line

ccf.2.line <- ggplot(dedup.maf.sbs.with.probable.process.ccf_hat |> filter(maxi_group=="2"),
                     aes(which_max_SBS, ccf_hat, fill=which_max_SBS)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale="width") +
  scale_fill_manual(values=sbs.palette) +
  facet_grid(~ which_max_SBS_group, scales = "free_x", space = "free_x") +
  labs(x="", y="Cancer Cell Fraction Maximum Estimate") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ccf.2.line

library(patchwork)
ccf.sbs.plot <- ccf.1.line + ccf.2.line + plot_layout(ncol = 1)
pdf("../figures/sbs.ccf.plot.over.50.ccfhat.pdf")
print(ccf.sbs.plot)
dev.off()

# avg per pt; clustered/nonclustered -------

average.dedup.sbs <- dedup.maf.sbs |>
  group_by(Participant_ID, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group) |>
  summarise_at(sbs.cols, sum) |>
  mutate(total= sum(c_across(any_of(sbs.cols)))) |>
  mutate_at(sbs.cols, ~./total) |>
  pivot_longer(-c(Participant_ID, total, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA")) |>
  mutate(SBS=factor(SBS, levels=rev(sbs.order))) |>
  mutate(Group=factor(Group, levels=rev(c("Genome", "clustered", "IG"))))

# average overall across patients
avg.average.dedup.sbs <- average.dedup.sbs |>
  group_by(Group, SBS) |>
  summarise(weight=sum(weight)/n())


# AVG SBS vs NEW ----------------------------------------------------------

avg.average.dedup.sbs |> 
  filter(Group=="Genome") |>
  mutate(COSMIC=str_detect(SBS, "SBS")) |>
  group_by(COSMIC) |> 
  summarise(weight=sum(weight))


# tim's plot --------------------------------------------------------------

ggplot(avg.average.dedup.sbs, aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )

### avg per pt; clonal; clustered/nonclustered -------

average.dedup.sbs.clonal <- dedup.maf.sbs |>
  group_by(Participant_ID, Gender, Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group, clonal_landau) |>
  summarise_at(sbs.cols, sum) |>
  mutate(total= sum(c_across(any_of(sbs.cols)))) |>
  mutate_at(sbs.cols, ~./total) |>
  pivot_longer(-c(Participant_ID, total, Gender, Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group, clonal_landau), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA")) |>
  mutate(SBS=factor(SBS, levels=rev(sbs.order))) |>
  mutate(Group=factor(Group, levels=rev(c("Genome", "clustered", "IG"))))

### duplicated clonal mutants -------

clonal.dup.hrd.mutations <- dedup.maf.sbs |>
  filter(clonal==1 & multiplicity==2 & Chromosome %in% c(3, 5, 7, 9, 11, 15, 19, 21)) |>
  group_by(Participant_ID, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group, clonal_landau) |>
  summarise_at(sbs.cols, sum) |>
  mutate(total= sum(c_across(any_of(sbs.cols)))) |>
  mutate_at(sbs.cols, ~./total) |>
  pivot_longer(-c(Participant_ID, total, Gender, Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group, clonal_landau), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA")) |>
  mutate(SBS=factor(SBS, levels=rev(sbs.order))) |>
  mutate(Group=factor(Group, levels=rev(c("Genome", "clustered", "IG"))))

avg.clonal.dup.hrd.mutations <- clonal.dup.hrd.mutations |>
  group_by(Group, SBS) |>
  summarise(weight=sum(weight)/n())
  # group_by(SBS, Group) |>
  # mutate(weight=weight / sum(weight))

### subclonal mutations -------

subclonal.dedup.sbs <- dedup.maf.sbs |>
  filter(non_clonal_landau==TRUE) |>
  group_by(Participant_ID,  Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group) |>
  summarise_at(sbs.cols, sum) |>
  mutate(total= sum(c_across(any_of(sbs.cols)))) |>
  mutate_at(sbs.cols, ~./total) |>
  pivot_longer(-c(Participant_ID,  Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, total, Group), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA")) |>
  mutate(SBS=factor(SBS, levels=rev(sbs.order))) |>
  mutate(Group=factor(Group, levels=rev(c("Genome", "clustered", "IG"))))

avg.subclonal.dedup.sbs <- subclonal.dedup.sbs |>
  group_by(Group, SBS) |>
  summarise(weight=sum(weight)/n())

## plot  -------------------------------------------------------

### show profile ------
ggplot(average.dedup.sbs |> filter(Participant_ID=="CSL141"), aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )

### show profile ------
ggplot(average.dedup.sbs.clonal |> filter(Participant_ID=="CSL141"), aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  ) +
  facet_wrap(~clonal_landau)

### overall sig sbs clonal duplicated in hrd chromosomes ----
overall.plot <- ggplot(avg.average.dedup.sbs, aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )

# overall sig sbs clonal duplicated in hrd chromosomes
hrd.dup.plot <- ggplot(avg.clonal.dup.hrd.mutations, aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )

# subclonal mutations
avg.subclonal.dedup.sbs.plot <- ggplot(avg.subclonal.dedup.sbs, aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )
# 
# ### subclonal cat ------
# avg.subclonal.cat.dedup.sbs.plot <- ggplot(avg.subclonal.cat.dedup.sbs |> filter(Group=="clustered" & !is.na(ccf_cat)), aes(x=weight, y=ccf_cat, fill=SBS)) +
#   geom_bar(stat="identity", position = "fill") +
#   scale_fill_manual(values=sbs.palette) +
#   theme_bw() +
#   theme(title = element_text(size=8),
#         axis.text.y = element_text(size=8),
#         axis.text.x = element_text(size=7),
#         panel.grid = element_blank(),
#         legend.position = "bottom"
#   )
# 
# avg.subclonal.cat.dedup.sbs.plot

sbs.plot.subclonal.clonal.dup <- overall.plot + theme(legend.position = "none") +
  hrd.dup.plot + theme(legend.position = "none") + avg.subclonal.dedup.sbs.plot +plot_layout(ncol = 1)

ggsave(filename = "../figures/sbs.plot.subclonal.clonal.dup.png", plot = sbs.plot.subclonal.clonal.dup, width = 7, height = 9)
ggsave(filename = "../figures/sbs.plot.subclonal.clonal.dup.pdf", plot = sbs.plot.subclonal.clonal.dup, width = 7, height = 9)



# CCF / clonality per mutation mechanism ----------------------------------

dedup.maf.sbs <- read_tsv("../data/_MM_Sigs_HDP/20240415_deduplicated_maf.tsv")

clonal.process.family <- dedup.maf.sbs |>
  filter(!startsWith(Tumor_Sample_Barcode, "MMRF")) |>
  filter(clonal_landau | non_clonal_landau) |>
  group_by(Guess, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Participant_ID) |>
  summarize(n=sum(clonal_landau), n_nonclonal=sum(non_clonal_landau), frac_clonal=n/(n+n_nonclonal))

# FIGURE 4B ---------------------------------------------------------------

clonal.process.family |>
  group_by(Guess) |>
  summarise(total_clonal=sum(n), total_mutations=sum(n)+sum(n_nonclonal), total_non_clonal=sum(n_nonclonal))
# # A tibble: 4 × 4
# Guess      total_clonal total_mutations total_non_clonal
# <chr>             <int>           <int>            <int>
# 1 AID/SHM           11415           15636             4221
# 2 APOBEC            49566          101922            52356
# 3 All<0.85         397876          733273           335397
# 4 Clock-like        14996           24780             9784

clonal.process.family |> 
  group_by(Guess) |> 
  summarize(median=median(frac_clonal))
# # A tibble: 4 × 2
# Guess      median
# <chr>       <dbl>
# 1 AID/SHM     0.714
# 2 APOBEC      0.3  
# 3 All<0.85    0.630
# 4 Clock-like  0.690

fracion_clonal.stats <- clonal.process.family |> 
  ungroup() |>
  filter(Guess!="All<0.85") |>
  mutate(Guess=factor(Guess, levels=c("Clock-like", "AID/SHM", "APOBEC"))) |>
  dunn_test(frac_clonal ~ Guess) |>
  add_xy_position() |>
  mutate(y.position=y.position+0.1)

# # A tibble: 3 × 13
# .y.         group1     group2     n1    n2 statistic        p    p.adj p.adj.signif y.position groups        xmin  xmax
# <chr>       <chr>      <chr>   <int> <int>     <dbl>    <dbl>    <dbl> <chr>             <dbl> <named list> <dbl> <dbl>
#   1 frac_clonal Clock-like AID/SHM   133   173      2.53 1.15e- 2 1.15e- 2 *                  1.10 <chr [2]>        1     2
# 2 frac_clonal Clock-like APOBEC    133    53     -5.28 1.32e- 7 2.64e- 7 ****               1.10 <chr [2]>        1     3
# 3 frac_clonal AID/SHM    APOBEC    173    53     -7.31 2.58e-13 7.75e-13 ****               1.10 <chr [2]>        2     3

signature.fraction.clonal.plot <- ggplot(clonal.process.family |> 
         filter(Guess!="All<0.85")
       , aes(factor(Guess, levels=c("Clock-like", "AID/SHM", "APOBEC"), labels=c("Clock-\nlike", "AID\n/SHM", "APO\nBEC")), 
             frac_clonal,
             color=factor(Guess, levels=c("Clock-like", "AID/SHM", "APOBEC")))) +
  geom_boxplot(size=0.25, outlier.size = 0.25) +
  # geom_jitter(alpha=0.2, color="black", shape=1) +
  labs(x="", y="Fraction of clonal mutations\n(Cancer cell fraction >85%)", fill="", color="") +
  scale_color_manual(values=c("Clock-like"="#B22222", "AID/SHM"="#1874CD", "APOBEC"="#FFA500")) +
  scale_y_continuous(labels=scales::percent, limits=c(0, NA), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  stat_pvalue_manual(fracion_clonal.stats, step.increase = 0.05, tip.length = 0, size = 2.5) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6))
signature.fraction.clonal.plot

ggsave(plot = signature.fraction.clonal.plot, 
       filename = "../figures/signature.fraction.clonal.plot.pdf", 
       width=1.4, height = 2, scale = 1)

# sbs.palette = 
c("N9"="#333333", "N12"="#5c5c5c", "N16"="#858585", "N18"="#adadad", "N19"="#d6d6d6", "N20"="#dddddd",
      "N13"="#C7E9C0", "N15"="#A1D99B", "SBS8"="#74C476", "SBS16"="#31A354", "SBS18"="#006D2C",
      "SBS1"="#FA8072", "SBS5"="#B22222", "SBS40"="#B22222", "SBS2"="#FFA500", "SBS13"="#EEEE00", 
      "SBS17a"="#EEAEEE", "SBS17b"="#B452CD", "SBS9"="#BFEFFF", "N17"="#00BFFF", "SBS84-like"="#00BFFF", "SBS85"="#1874CD")
  
  
## clonal.process.family.per.group -------
clonal.process.family.per.group <- dedup.maf.sbs |>
  filter(!startsWith(Tumor_Sample_Barcode, "MMRF")) |>
  filter(clonal_landau | non_clonal_landau) |>
  group_by(Guess, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Participant_ID, Group) |>
  summarize(n=sum(clonal_landau), frac_clonal=n/(n+sum(non_clonal_landau)))

ggplot(clonal.process.family.per.group, aes(Guess, frac_clonal, fill=Guess)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~Group) +
  stat_compare_means()


## clonal.process.family.per.IMWG -------

clonal.process.family.per.group <- dedup.maf.sbs |>
  filter(!startsWith(Tumor_Sample_Barcode, "MMRF")) |>
  filter(clonal_landau | non_clonal_landau) |>
  group_by(Guess, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Participant_ID) |>
  summarize(n=sum(clonal_landau), frac_clonal=n/(n+sum(non_clonal_landau)))

ggplot(clonal.process.family, aes(SimpleStage, frac_clonal, fill=IMWG)) +
  geom_boxplot() +
  theme_bw() +
  stat_compare_means() +
  facet_wrap(~Guess)

ggplot(clonal.process.family, aes(Guess, frac_clonal, fill=IMWG)) +
  geom_boxplot() +
  theme_bw() +
  stat_compare_means()



# dup Apr 23 -----------------------------------------------------------------


dup.clonal <- read_tsv("../data/_MM_Sigs_HDP/maf_patient_dup_all_2024_04_16.txt")

dedup.dup.clonal.maf.sbs <- tibble(dup.clonal) |>
  # filter(Participant_ID=="CTF001") |>
  group_by(Participant_ID, Chromosome, Start_position) |>
  slice_head(n=1) |>
  ungroup() |>
  mutate(Group = case_when(
    str_detect(SampleID, "_clustered$") ~ "clustered",
    str_detect(SampleID, "_IG$") ~ "IG",
    TRUE ~ "Genome")) |>
  inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID")) |>
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) |>
  mutate(SimpleStage=factor(Stage, levels=c("MGUS", "SMM", "MM"), labels=c("MGUS/SMM", "MGUS/SMM", "MM"))) |>
  mutate(Guess=case_when(
    SBS9 + SBS85 + N17 > 0.85  ~ "AID/SHM",
    SBS2 + SBS13 > 0.85 ~ "APOBEC",
    SBS1 + SBS5 + SBS40 > 0.85 ~ "Clock-like",
    TRUE ~ "All<0.85")) |>
  mutate(clonal_landau = ccf_median>=0.85 & !is.na(ccf_median), non_clonal_landau= ccf_median<0.85 & !is.na(ccf_median))


average.dedup.sbs <- dedup.dup.clonal.maf.sbs |>
  group_by(Participant_ID, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group) |>
  summarise_at(sbs.cols, sum) |>
  mutate(total= sum(c_across(any_of(sbs.cols)))) |>
  mutate_at(sbs.cols, ~./total) |>
  pivot_longer(-c(Participant_ID, total, Gender,Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020, HDP_Reference_Pair, SimpleStage, Group), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=fct_recode(SBS, `SBS84-like`="N17")) |>
  mutate(SBS.group=case_when(
    SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20") ~ "Possible sequencing artefacts",  
    SBS %in% c("N13", "N15") ~ "Novel, Unknown",  
    SBS %in% c("SBS8", "SBS16") ~ "Unknown",
    SBS %in% c("SBS18") ~ "ROS",
    SBS %in% c("SBS1", "SBS5", "SBS40") ~ "Clock-like",
    SBS %in% c("SBS2", "SBS13") ~ "APOBEC",
    SBS %in% c("SBS17a","SBS17b") ~ "Unknown, 17",
    SBS %in% c("SBS9", "SBS84-like","SBS85") ~ "AID/SHM",
    TRUE ~ "BLABLA")) |>
  mutate(SBS=factor(SBS, levels=rev(sbs.order))) |>
  mutate(Group="Clonal duplicated before Hyperdiploidy")

avg.clonal.dup.sbs <- average.dedup.sbs |>
  group_by(Group, SBS) |>
  summarise(weight=sum(weight)/n())

avg.clonal.dup.sbs.plot <- ggplot(avg.clonal.dup.sbs, aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )


# clonal dup, genome, clustered, ig mutations -----------------------------

overall.clonal.dup.genome.clustered.ig.plot <- bind_rows(avg.average.dedup.sbs, avg.clonal.dup.sbs) |>
  mutate(Group=factor(Group, 
                      levels=rev(c("Clonal duplicated before Hyperdiploidy", "Genome", "clustered", "IG")),
                      labels=rev(c("Early mutations", "Genome-wide\nlandscape", "Clustered\nmutations", "Immunoglobulin\nmutations")))) |>
  ggplot(aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill", color="#333333", linewidth=0.25) +
  scale_fill_manual(values=sbs.palette) +
  theme_bw() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        legend.position = "bottom"
  )

ggsave("../figures/avg.clonal.dup.sbs.png", avg.clonal.dup.sbs.plot, width = 6, height = 4)
ggsave("../figures/avg.clonal.dup.sbs.pdf", avg.clonal.dup.sbs.plot, width = 6, height = 4)

ggsave("../figures/overall.clonal.dup.genome.clustered.ig.plot.png", overall.clonal.dup.genome.clustered.ig.plot, width = 6, height = 4)
ggsave("../figures/overall.clonal.dup.genome.clustered.ig.plot.pdf", overall.clonal.dup.genome.clustered.ig.plot, width = 6, height = 4)


# May 24 24 ---------------------------------------------------------------

# Extract only duplicated and compare with genome wide
ALL <- bind_rows(avg.average.dedup.sbs, avg.clonal.dup.sbs) |>
  mutate(Group=factor(Group, 
                      levels=rev(c("Clonal duplicated before Hyperdiploidy", "Genome", "clustered", "IG")),
                      labels=rev(c("Early mutations", "Genome-wide\nlandscape", "Clustered\nmutations", "Immunoglobulin\nmutations"))))

write_tsv(ALL, "../data/_MM_Sigs_HDP/20240524_export_weights_with_duplicated.tsv")



overall.clonal.dup.genome.plot <- bind_rows(avg.average.dedup.sbs, avg.clonal.dup.sbs) |>
  mutate(Group=factor(Group, 
                      levels=rev(c("Clonal duplicated before Hyperdiploidy", "Genome", "clustered", "IG")),
                      labels=rev(c("Early mutations", "Genome-wide\nlandscape", "Clustered\nmutations", "Immunoglobulin\nmutations")))) |>
  filter(Group %in% c("Early mutations", "Genome-wide\nlandscape")) |>
  ggplot(aes(x=weight, y=Group, fill=SBS)) +
  geom_bar(stat="identity", position = "fill", color="#333333", linewidth=0.25) +
  scale_fill_manual(values=sbs.palette) +
  theme_nothing() +
  theme(title = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=7),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none"
  )


bind_rows(avg.average.dedup.sbs, avg.clonal.dup.sbs) |>
  mutate(Group=factor(Group, 
                      levels=rev(c("Clonal duplicated before Hyperdiploidy", "Genome", "clustered", "IG")),
                      labels=rev(c("Early mutations", "Genome-wide\nlandscape", "Clustered\nmutations", "Immunoglobulin\nmutations")))) |>
  filter(Group %in% c("Early mutations", "Genome-wide\nlandscape")) |>
  chisq_test(x=)

overall.clonal.dup.genome.plot

ggsave("../figures/avg.clonal.dup.genome.sbs.png", overall.clonal.dup.genome.plot, width = 6, height = 4)
ggsave("../figures/avg.clonal.dup.genome.sbs.pdf", overall.clonal.dup.genome.plot, width = 6, height = 4)



avg.average.dedup.sbs |> 
  filter(Group=="Genome") |> 
  filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20"))

avg.average.dedup.sbs |> 
  filter(Group=="Genome") |> 
  filter(SBS %in% c("N9", "N12", "N16", "N18", "N19", "N20")) |>
  summarise(sum(weight))
