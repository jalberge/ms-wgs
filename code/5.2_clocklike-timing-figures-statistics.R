setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(readxl)

source("0_annotate_samples.R")

initial.disease.stage.palette <- c("MM"="#253494", 
                                   "MGUS/SMM\nprogressors"="#225ea8", 
                                   "SMM"="#41b6c4", 
                                   "MGUS"="#7fcdbb")

# better
disease.progression.palette <- c("MGUS\nStable"="#9AB6C6", 
                                 "MGUS\nProgressed"="#316B9B", 
                                 "SMM\nStable"="#A3C487", 
                                 "SMM\nProgressed"="#338E2E", 
                                 "MM"="#D58E8D")
disease.stage.palette <- c("MGUS"="#9AB6C6", 
                           "SMM"="#A3C487", 
                           "MM"="#D58E8D")


# tsv ---------------------------------------------------------------------

maura <- read_tsv('../data/Maura_data_MM_timing_2024_03_13.txt')
maura.annot <- read_xlsx('../data/Maura_data_MM_clinical_data.xlsx')

head(maura)
head(maura.annot)

table(maura$Sample %in% maura.annot$Samples)
maura$Sample[!(maura$Sample %in% maura.annot$Samples)]

# viz ---------------------------------------------------------------------

maura.timing <- maura |> 
  filter(!is.na(as.numeric(HRD_timing))) |>
  mutate(across(HRD_timing, as.numeric)) |>
  inner_join(maura.annot, by=c("Sample"="Samples")) |>
  mutate(Age=as.numeric(`Age (years)`)) |>
  mutate(
    Absolute_HRD_timing=(HRD_timing)*Age,
    Absolute_HRD_timing_low=(HRD_timing_low)*Age,
    Absolute_HRD_timing_high=(HRD_timing_high)*Age
  ) |>
  group_by(Patients) |> 
  mutate(Timepoint=paste0("T", row_number())) |>
  mutate(Initial_stage=Stage[1]) |>
  mutate(Initial_stage=factor(Initial_stage, levels=c("SMM - progressed", "MM"))) |>
  mutate(Patients = fct_reorder(as.factor(Patients), Age))


# Relative scale ----------------------------------------------------------

ggplot(maura.timing, aes(x=Patients, y=HRD_timing, group=Sample, color=Stage)) +
  geom_point() +
  geom_linerange(aes(ymin=HRD_timing_low, ymax=HRD_timing_high))

# Absolute scale ----------------------------------------------------------

# factor to shift different samples from the same patient on the x axis
dodge.factor <- 0.7

absolute.age <- ggplot(maura.timing, aes(x=fct_reorder(Patients, Age, min), y=Absolute_HRD_timing, color=Stage, group=Sample)) +
  # geom_segments does not allow dodging; however geom_segment does :-)  https://stackoverflow.com/a/21922792/4783389
  # geom_segment(aes(x=Patients, xend=Patients, y=Absolute_HRD_timing, yend=Age), color="grey", linetype=2, position=position_dodge(width = dodge.factor)) +
  ### data
  geom_linerange(aes(ymin=Absolute_HRD_timing, ymax=Age), color="grey", linetype=2, position=position_dodge(width = dodge.factor)) +
  geom_point(aes(y=Age), color="black", position=position_dodge(width = dodge.factor)) +
  geom_point(position=position_dodge(width = dodge.factor), shape=21) +
  geom_linerange(aes(ymin=Absolute_HRD_timing_low, ymax=Absolute_HRD_timing_high), position=position_dodge(width = dodge.factor)) +
  ### scales
  scale_color_manual(values=c("MM"="darkblue", "MGUS/SMM\nprogressors"="royalblue")) +
  scale_y_continuous(limits = c(0, NA)) +
  ### organization
  facet_grid(~Initial_stage, scales = "free_x", space="free_x") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        legend.position = "top", strip.background = element_blank(), ) +
  ### annot
  labs(x="MGUS/SMM/MM patients", y="Age (years)", color="Disease stage")
absolute.age

ggsave(plot = absolute.age, filename = "../figures/timing_maura_samples.png", width = 4, height = 2, scale = 1.5)

# su2c ---------------------------------------------------------

local.t <- read_tsv("../data/MM_HRD_timing_all_2024_03_06.txt")

local.timing <- local.t |>
  filter(!is.na(as.numeric(HRD_timing))) |>
  mutate(across(HRD_timing, as.numeric)) |>
  inner_join(clinical.and.terra.ref.pairs.samples, by=c("Sample ID"="HDP_Participant_ID")) |>
  mutate(Age=as.numeric(Age.x)) |>
  mutate(
    Absolute_HRD_timing=(HRD_timing)*Age,
    Absolute_HRD_timing_low=(HRD_timing_low)*Age,
    Absolute_HRD_timing_high=(HRD_timing_high)*Age
  ) |>
  mutate(Stage=case_when(
    Stage %in% c("MGUS")~"MGUS",
    Stage %in% c("LRSMM", "IRSMM", "HRSMM")~"SMM",
    Stage %in% c("MM", "NDMM")~"MM")
                         ) |>
  mutate(Initial_stage=Stage) |>
  mutate(Sample=`Sample ID`) |>
  group_by(`Sample ID`) |> 
  mutate(Timepoint=paste0("T", row_number())) |>
  mutate(Patients = fct_reorder(as.factor(`Sample ID`), Age))


## Relative scale ----------------------------------------------------------

ggplot(local.timing, aes(x=Patients, y=HRD_timing, group=Sample, color=Stage)) +
  geom_point() +
  geom_linerange(aes(ymin=HRD_timing_low, ymax=HRD_timing_high))

## absolute age

# factor to shift different samples from the same patient on the x axis
dodge.factor <- 0.7

absolute.local.age <- ggplot(local.timing, aes(x=fct_reorder(Patients, Age, min), y=Absolute_HRD_timing, color=Stage, group=Sample)) +
  ### data
  geom_linerange(aes(ymin=Absolute_HRD_timing, ymax=Age), color="grey", linetype=2, position=position_dodge(width = dodge.factor)) +
  geom_point(aes(y=Age), color="black", position=position_dodge(width = dodge.factor)) +
  geom_point(position=position_dodge(width = dodge.factor), shape=21) +
  geom_linerange(aes(ymin=Absolute_HRD_timing_low, ymax=Absolute_HRD_timing_high), position=position_dodge(width = dodge.factor)) +
  ### scales
  scale_color_manual(values=c("MM"="darkblue", "SMM"="lightblue", "MGUS"="lightgreen")) +
  scale_y_continuous(limits = c(0, NA)) +
  ### organization
  facet_grid(~Initial_stage, scales = "free_x", space="free_x") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        legend.position = "top", strip.background = element_blank(), ) +
  ### annot
  labs(x="MGUS/SMM/MM patients", y="Age (years)", color="Disease stage")
absolute.local.age

ggsave(plot = absolute.local.age, filename = "../figures/timing_local.png", width = 4, height = 2, scale = 1.5)


# Maura and SU2C ----------------------------------------------------

all.timing <- bind_rows(maura.timing, local.timing, .id = "source")

# unfortunately ids are from redcap and need to remain controlled-access. 
# TODO: switch to publication ID (SCIxxx or CTFxxx) 
source("5.0-progressor-ids.R")

all.timing.groups <- all.timing |>
  group_by(Patients) |>
  slice_head(n=1) |>

  mutate(Disease_stage = ifelse(Stage=="SMM - progressed", "SMM", Stage)) |>
  mutate(Simple_Disease_stage = ifelse(Disease_stage %in% c("MGUS", "SMM"), "MGUS/SMM", Disease_stage)) |>
  mutate(Progression_status = ifelse(Participant_ID %in% c(su2c.progressors.w.timing, mgus.oben.prog) | Stage=="SMM - progressed", "Progressed", "Stable")) |>
  mutate(Disease_Progression = ifelse(Disease_stage %in% c("MGUS", "SMM"), paste0(Disease_stage, "\n", Progression_status), Disease_stage)) |>
  mutate(Simple_Disease_Progression = ifelse(Simple_Disease_stage == "MGUS/SMM", paste0(Simple_Disease_stage, "\n", Progression_status), Simple_Disease_stage)) |>
  
  mutate(Disease_stage = factor(Disease_stage, levels=c("MGUS", "SMM", "MM"))) |>
  mutate(Simple_Disease_stage = factor(Simple_Disease_stage, levels=c("MGUS/SMM", "MM"))) |>
  mutate(Progression_status = factor(Progression_status, levels=c("Stable", "Progressed"))) |>
  mutate(Disease_Progression = factor(Disease_Progression, levels=c("MGUS\nStable", "MGUS\nProgressed", "SMM\nStable", "SMM\nProgressed", "MM"))) |>
  mutate(Simple_Disease_Progression = factor(Simple_Disease_Progression, levels=c("MGUS/SMM\nStable", "MGUS/SMM\nProgressed", "MM"))) |>
  # mutate(Initial_stage_fct = fct_recode(Initial_stage_fct, "MGUS/SMM\nprogressors"="SMM - progressed")) |>
  
  mutate(
    Tumor_Age=Age-Absolute_HRD_timing, .after = purity,
    Tumor_Age_high=Age-Absolute_HRD_timing_low,
    Tumor_Age_low=Age-Absolute_HRD_timing_high,
         ) |>
  ungroup() |>
  select(source:Initial_stage, 
         Disease_stage, Progression_status, Disease_Progression, Simple_Disease_stage, Simple_Disease_Progression)

# SAVE  ---------

all.timing.groups |> write_tsv("../data/Timing_MGUS_SMM_MGUS_Progression_20230401.tsv", quote = "needed")

# AGE OF ONSET -----------------------------------------------------------

## Detailed ----------------

absolute.age <- ggplot(all.timing.groups, aes(x=fct_reorder(Patients, Age, min), y=Absolute_HRD_timing, color=Disease_Progression, group=Sample)) +
  ### data
  geom_linerange(aes(ymin=Absolute_HRD_timing, ymax=Age), color="grey", linetype=2, position=position_dodge(width = dodge.factor)) +
  geom_point(aes(y=Age), color="black", position=position_dodge(width = dodge.factor)) +
  geom_point(position=position_dodge(width = dodge.factor), shape=21) +
  geom_linerange(aes(ymin=Absolute_HRD_timing_low, ymax=Absolute_HRD_timing_high), position=position_dodge(width = dodge.factor)) +
  ### scales
  scale_color_manual(values=darken(brewer.pal(6, "Paired"))) +
  scale_y_continuous(limits = c(0, NA)) +
  ### organization
  facet_grid(~Disease_stage, scales = "free_x", space="free_x") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.background = element_blank(), ) +
  ### annot
  labs(x="", y="Age (years)", color="Disease stage")

absolute.age

ggsave(plot = absolute.age, filename = "../figures/timing_total.png", width = 4, height = 2, scale = 1.5)
ggsave(plot = absolute.age, filename = "../figures/timing_total.pdf", width = 4, height = 2, scale = 1.5)


## Boxplot MGUS/SMM/MM ------------------

all.timing.groups |> kruskal_test(Absolute_HRD_timing~Disease_stage)
all.timing.groups |> dunn_test(Absolute_HRD_timing~Disease_stage)
all.timing.groups |> group_by(Disease_stage) |> summarise(median=median(Absolute_HRD_timing))

p.vals.stage <- all.timing.groups |> dunn_test(Absolute_HRD_timing~Disease_stage) |> add_xy_position()

boxplot.age.stage <- ggplot(all.timing.groups, aes(Disease_stage, Absolute_HRD_timing, color=Disease_stage)) +
  geom_boxplot() +
  geom_jitter(alpha=0.3) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank()) +
  scale_color_manual(values=disease.stage.palette) +
  labs(x="", y="Estimated age of onset of the disease (years)", color="Disease stage") +
  stat_pvalue_manual(p.vals.stage, hide.ns = TRUE, tip.length = 0) +
  stat_compare_means(label.x.npc = "left", label.y.npc = "bottom", size=2.5)

boxplot.age.stage

ggsave(plot = boxplot.age.stage, filename = "../figures/boxplot_timing_age.png", width = 1.5, height = 2, scale = 1.5)
ggsave(plot = boxplot.age.stage, filename = "../figures/boxplot_timing_age.pdf", width = 4, height = 2, scale = 1.5)


## Boxplot Stable/prog ------------------

mgus.smm.timing.groups <- all.timing.groups |> filter(Simple_Disease_stage!="MM")

mgus.smm.timing.groups |> wilcox_test(Absolute_HRD_timing~Progression_status)
p.vals.progression <- mgus.smm.timing.groups |> dunn_test(Absolute_HRD_timing~Progression_status) |> add_xy_position()
mgus.smm.timing.groups |> group_by(Progression_status) |> summarise(median=median(Absolute_HRD_timing))

boxplot.age.progression <- ggplot(mgus.smm.timing.groups, aes(Progression_status, Absolute_HRD_timing, color=Progression_status)) +
  geom_boxplot() +
  geom_jitter(alpha=0.3) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank()) +
  scale_color_manual(values=brewer.pal(3, "Set1")) +
  labs(x="", y="Estimated age of onset of the disease (years)", color="Disease stage") +
  stat_compare_means(label = "p", label.x.npc = "right", label.y.npc = "bottom", size=2.5, hide.ns = TRUE, tip.length = 0) +
  stat_pvalue_manual(p.vals.progression, hide.ns = TRUE, tip.length = 0, label = "p.adj.signif")
  
boxplot.age.progression

ggsave(plot = boxplot.age.progression, filename = "../figures/boxplot_timing_age_progression.png", width = 1.5, height = 2, scale = 1.5)
ggsave(plot = boxplot.age.progression, filename = "../figures/boxplot_timing_age_progression.pdf", width = 4, height = 2, scale = 1.5)


## Group -------------------------------------------------------------------

timing.details.plus.average <- plot_grid(absolute.age, boxplot.age.stage, boxplot.age.progression, rel_widths = c(3,1,0.8), align = "h", axis = "tb", nrow = 1)

timing.details.plus.average

ggsave("../figures/timing.details.plus.boxplots.png", timing.details.plus.average, width = 6, height = 2.4, scale = 1.5)
ggsave("../figures/timing.details.plus.boxplots.pdf", timing.details.plus.average, width = 6, height = 2.4, scale = 1.5)


## stats -----

dim(all.timing.groups)
table(all.timing.groups$Disease_stage)





# TUMOR AGE -----------------------------------------------------------

## Detailed ----------------

absolute.tumor.age <- ggplot(all.timing.groups, aes(x=fct_reorder(Patients, Tumor_Age, min), y=Tumor_Age, color=Disease_Progression, group=Sample)) +
  geom_point(position=position_dodge(width = dodge.factor), shape=21) +
  geom_linerange(aes(ymin=Tumor_Age_low, ymax=Tumor_Age_high), position=position_dodge(width = dodge.factor)) +
  ### scales
  scale_color_manual(values=darken(brewer.pal(6, "Paired"))) +
  scale_y_continuous(limits = c(NA, NA)) +
  ### organization
  facet_grid(~Disease_stage, scales = "free_x", space="free_x") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), 
        # axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.background = element_blank(), ) +
  ### annot
  labs(x="", y="Tumor Age (years)", color="Disease stage")

absolute.tumor.age

ggsave(plot = absolute.tumor.age, filename = "../figures/supp_timing_total.png", width = 4, height = 2, scale = 1.5)
ggsave(plot = absolute.tumor.age, filename = "../figures/supp_timing_total.pdf", width = 4, height = 2, scale = 1.5)


## Boxplot MGUS/SMM/MM ------------------

all.timing.groups |> kruskal_test(Tumor_Age~Disease_stage)
# # A tibble: 1 × 6
# .y.           n statistic    df       p method        
# * <chr>     <int>     <dbl> <int>   <dbl> <chr>         
#   1 Tumor_Age    91      11.0     2 0.00417 Kruskal-Wallis
all.timing.groups |> dunn_test(Tumor_Age~Disease_stage)
# .y.       group1 group2    n1    n2 statistic       p   p.adj p.adj.signif
# * <chr>     <chr>  <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
#   1 Tumor_Age MGUS   SMM       14    46      1.56 0.120   0.120   ns          
# 2 Tumor_Age MGUS   MM        14    31      3.14 0.00170 0.00511 **          
#   3 Tumor_Age SMM    MM        46    31      2.31 0.0211  0.0423  *   
all.timing.groups |> group_by(Disease_stage) |> summarise(median=median(Tumor_Age))
# # A tibble: 3 × 2
# Disease_stage median
# <fct>          <dbl>
#   1 MGUS            16.2
# 2 SMM             20.5
# 3 MM              32.8
p.vals.stage <- all.timing.groups |> dunn_test(Tumor_Age~Disease_stage) |> add_xy_position()

boxplot.age.stage.tumor.age <- ggplot(all.timing.groups, aes(Disease_stage, Tumor_Age, color=Disease_stage)) +
  geom_boxplot() +
  geom_jitter(alpha=0.3) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank()) +
  scale_color_manual(values=disease.stage.palette) +
  labs(x="", y="Estimated tumor age (years)", color="Disease stage") +
  stat_pvalue_manual(p.vals.stage, hide.ns = TRUE, tip.length = 0) +
  stat_compare_means(label.x.npc = "middle", label.y.npc = "bottom", size=2.5)

boxplot.age.stage.tumor.age

ggsave(plot = boxplot.age.stage.tumor.age, filename = "../figures/supp_boxplot_timing_age.png", width = 1.5, height = 2, scale = 1.5)
ggsave(plot = boxplot.age.stage.tumor.age, filename = "../figures/boxplot_timing_age.pdf", width = 4, height = 2, scale = 1.5)


## Boxplot Stable/prog ------------------

mgus.smm.timing.groups <- all.timing.groups |> filter(Simple_Disease_stage!="MM")

#### REVISIONS #1 ####
# revisions: is age at HRD an independent factor of progression
independent.hrd.prog <- mgus.smm.timing.groups |> 
  inner_join(clinical.and.terra.ref.pairs.samples, by=c("Sample"="HDP_Participant_ID"))
library(logistf) # for Firth correction
# run 1 model without MGUS which is "cleaner" from a clinical point of view (we could get criticized that MGUS will not progress anyway)
independent.hrd.prog$Disease_Status <- factor(independent.hrd.prog$Disease_Status, levels=c("MGUS", "LRSMM", "IRSMM", "HRSMM"))
independent.hrd.prog$Tumor_Age_Decade <- independent.hrd.prog$Tumor_Age/10
# model:
hrd.age.independent.model <- glm(Progression_status~Tumor_Age_Decade + Disease_Status, data=independent.hrd.prog |> filter(Disease_Status!="MGUS"), family = "binomial")
#coefficients
tidy_model <- tidy(hrd.age.independent.model, conf.int = TRUE) |>
  mutate(estimate = exp(estimate), # exponentiate to have odds ratio scale
         conf.low = exp(conf.low),
         conf.high = exp(conf.high)) |>
  mutate(# trim at 5 to make the graph readable
         conf.high=pmin(conf.high, 5))
tidy_model
# A tibble: 4 × 7
# term                estimate std.error statistic p.value conf.low conf.high
# <chr>                  <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 (Intercept)           0.0262     1.36    -2.67   0.00758 0.000866     0.245
# 2 Tumor_Age_Decade      1.95       0.329    2.04   0.0415  1.08         4.11 
# 3 Disease_StatusIRSMM   0.940      1.43    -0.0429 0.966   0.0530       5    
# 4 Disease_StatusHRSMM   0.946      1.35    -0.0407 0.968   0.0692       5    

#forest
hrd.plot <- ggplot(tidy_model |> filter(term!="(Intercept)"), aes(x = estimate, y = term, color=p.value<0.05)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  theme_minimal() +
  xlab("Odds Ratio (Exponentiated Coefficients)") +
  ylab("") +
  scale_color_manual(values=c("TRUE"=brewer.pal(4, "Set1")[2], "FALSE"="#777777")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#333333") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
ggsave("../figures/hrd.age.independent.pdf", width = 3, height = 3)
# 
# term                estimate std.error statistic p.value conf.low conf.high
# <chr>                  <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#   1 (Intercept)           0.0262     1.36    -2.67   0.00758 0.000866     0.245
# 2 Tumor_Age_Decade      1.95       0.329    2.04   0.0415  1.08         4.11 
# 3 Disease_StatusIRSMM   0.940      1.43    -0.0429 0.966   0.0530      25.4  
# 4 Disease_StatusHRSMM   0.946      1.35    -0.0407 0.968   0.0692      23.5  

#### REVISIONS ####


mgus.smm.timing.groups |> wilcox_test(Tumor_Age~Progression_status)
# # A tibble: 1 × 7
# .y.       group1 group2        n1    n2 statistic       p
# * <chr>     <chr>  <chr>      <int> <int>     <dbl>   <dbl>
#   1 Tumor_Age Stable Progressed    44    16       192 0.00677
p.vals.progression <- mgus.smm.timing.groups |> dunn_test(Tumor_Age~Progression_status) |> add_xy_position()
# A tibble: 1 × 13
# .y.       group1 group2        n1    n2 statistic       p   p.adj p.adj.signif y.position groups        xmin  xmax
# <chr>     <chr>  <chr>      <int> <int>     <dbl>   <dbl>   <dbl> <chr>             <dbl> <named list> <dbl> <dbl>
#   1 Tumor_Age Stable Progressed    44    16      2.67 0.00748 0.00748 **                 61.4 <chr [2]>        1     2
mgus.smm.timing.groups |> group_by(Progression_status) |> summarise(median=median(Tumor_Age))
# # A tibble: 2 × 2
# Progression_status median
# <fct>               <dbl>
#   1 Stable               15.4
# 2 Progressed           29.6
boxplot.age.progression.tumor.age <- ggplot(mgus.smm.timing.groups, aes(Progression_status, Tumor_Age, color=Progression_status)) +
  geom_boxplot() +
  geom_jitter(alpha=0.3) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", panel.grid = element_blank()) +
  scale_color_manual(values=brewer.pal(3, "Set1")) +
  labs(x="", y="Estimated tumor age (years)", color="Disease stage") +
  stat_compare_means(label = "p", label.x.npc = "right", label.y.npc = "bottom", size=2.5, hide.ns = TRUE, tip.length = 0) +
  stat_pvalue_manual(p.vals.progression, hide.ns = TRUE, tip.length = 0, label = "p.adj.signif")

boxplot.age.progression.tumor.age

ggsave(plot = boxplot.age.progression.tumor.age, filename = "../figures/supp_boxplot_timing_age_progression.png", width = 1.5, height = 2, scale = 1.5)
ggsave(plot = boxplot.age.progression.tumor.age, filename = "../figures/supp_boxplot_timing_age_progression.pdf", width = 4, height = 2, scale = 1.5)


## Group -------------------------------------------------------------------

timing.details.plus.average.tumor.age <- plot_grid(absolute.tumor.age, boxplot.age.stage.tumor.age, boxplot.age.progression.tumor.age, rel_widths = c(3,1,0.8), align = "h", axis = "tb", nrow = 1)

timing.details.plus.average.tumor.age

ggsave("../figures/supp_timing.details.plus.boxplots.png", timing.details.plus.average.tumor.age, width = 6, height = 2.4, scale = 1.5)
ggsave("../figures/supp_timing.details.plus.boxplots.pdf", timing.details.plus.average.tumor.age, width = 6, height = 2.4, scale = 1.5)

## stats -----

dim(all.timing.groups)
table(all.timing.groups$Disease_stage)

# Mix age of onset / tumor age

timing.details.plus.average.tumor.age <- plot_grid(absolute.tumor.age, boxplot.age.stage.tumor.age, boxplot.age.progression.tumor.age, rel_widths = c(3,1,0.8), align = "h", axis = "tb", nrow = 1)

timing.details.plus.average.tumor.age

ggsave("../figures/supp_timing.details.plus.boxplots.png", timing.details.plus.average.tumor.age, width = 6, height = 2.4, scale = 1.5)
ggsave("../figures/supp_timing.details.plus.boxplots.pdf", timing.details.plus.average.tumor.age, width = 6, height = 2.4, scale = 1.5)
