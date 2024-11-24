setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(survival)
library(viridis)
library(survminer)

source("0_annotate_samples.R")

# NOTES FROM 2024-03-26
# RE-calculating score to match final matrix from Fig 1A, and add more clinical data


# Compute score again -----------------------------------------------------

# this is an important part of the paper so I'll do everything from scratch again with the final matrix.

list.files(path = "../data/", pattern = "FINAL_.*")

# matrix
# !!
# still some RedCap IDs in these sheets - can't make public on github
mm.meta <- read_tsv("../data/FINAL_meta_columns_matrix.tsv")

mm.t <- read_tsv("../data/FINAL_Jan4_TX_Focal_Numeric_Matrix.tsv")
mm.cnv <- read_tsv("../data/FINAL_Nov2_CNV_Focal_Numeric_Matrix.tsv")
mm.hrd <- read_tsv("../data/FINAL_Nov2_HRD_Counts_Numeric_Matrix.tsv")
mm.snv <- read_tsv("../data/FINAL_Oct9_SNV_Hotspot_Numeric_Matrix.tsv")
mm.sv <- read_tsv("../data/FINAL_Oct9_SV_Hotspot_Numeric_Matrix.tsv")

# these ones are not in the final analysis
mm.meta <- mm.meta |> filter(!(`Sample ID`%in%c("MMRF_1889", "MMRF_2363")))


mm.t.1 <- mm.t %>% 
  # select(Variable, starts_with("MMRF")) %>%
  pivot_longer(cols = !Variable) %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  rename(ID=name,
         t14_16="t(14;16)(MAF)",
         t14_20="t(14;20)(MAFB)",
         MYC_translocation="t(MYC;IgH/K/L)")

mm.cnv.1 <- mm.cnv %>% 
  # select(Variable, starts_with("MMRF")) %>%
  pivot_longer(cols = !Variable) %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  rename(ID=name,
         amp_G.peak.4_amp_16p13.13="+16p13.13",
         amp_chr_1q="+1q",
         amp_G.peak.1_amp_1q21.2="+1q21.2",
         amp_chr_6p="+6p", 
         del_G.peak.6_del_1p22.1="-1p22.1",
         del_G.peak.7_del_1p12="-1p12",
         del_G.peak.5_del_1p32.3="-1p32.3",
         del_G.peak.30_del_16q12.1="-16q12.1",
         del_chr_16q="-16q",
         del_G.peak.11_del_4p16.3="-4p16.3",
         del_chr_8p="-8p", 
         del_G.peak.17_del_8q24.21="-8q24.21a",
         del_G.peak.18_del_8q24.21="-8q24.21b",
         del_G.peak.16_del_8p23.3="-8p23.3",
         del_G.peak.27_del_14q23.3="-14q23.3",
         del_G.peak.28_del_14q31.1="-14q31.1", 
         del_G.peak.31_del_17q21.2="-17q21.2",
         amp_G.peak.2_amp_3q26.2="+3q26.2", 
         amp_G.peak.3_amp_8q24.21="+8q24.21", 
         del_G.peak.22_del_12p13.33="-12p13.33",
         del_G.peak.12_del_4q34.3="-4q34.3")

mm.hrd.1 <- mm.hrd %>% 
  # select(Variable, starts_with("MMRF")) %>%
  pivot_longer(cols = !Variable) %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  rename(ID=name,
         HyperDiploidy=HRD_chr_count) %>%
  mutate(HyperDiploidy=if_else(HyperDiploidy >=2, 1,0))

mm.snv.1 <- mm.snv %>% 
  # select(Variable, starts_with("MMRF")) %>%
  pivot_longer(cols = !Variable) %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  rename(ID=name)

colnames(mm.snv.1)[c(2:47)] <- paste0("SNV_", colnames(mm.snv.1)[c(2:47)])


mm <- mm.meta |>
  left_join( mm.hrd.1, by=c("Sample ID"="ID")) |>
  left_join( mm.snv.1, by=c("Sample ID"="ID")) |>
  left_join( mm.t.1, by=c("Sample ID"="ID")) |>
  left_join( mm.cnv.1, by=c("Sample ID"="ID"))
# mm <- left_join(mm.cnv.1, mm.hrd.1, by="ID")
# mm <- left_join(mm, mm.snv.1, by="ID")
# mm <- left_join(mm, mm.t.1, by=c("ID" = "Sample ID")

write_tsv(mm, "../data/export_MM-like_features_MGUS_SMM_MM.tsv")

# Specific for TP53 reprocessing ------------------------------------------

TP53 <- read_tsv("../JB/Original data/CN_TP53_all_datasets.txt")

TP53$ID <- gsub("_1_BM_CD138pos_pair", "", TP53$ID)
TP53$ID <- str_replace_all(TP53$ID, "-", "_")
TP53$ID <- str_replace_all(TP53$ID, "_pair", "_Tumor")
# SMM-064_pair
# SMM_064_Tumor
andrea.samples <- sort(unique(TP53$ID[!startsWith(TP53$ID, "MMRF") & !startsWith(TP53$ID, "Tumor_")& !startsWith(TP53$ID, "SMM_")]))
# > andrea.samples[andrea.samples %nin% clinical.and.terra.ref.pairs.samples$Reference_Pair]
# [1] "CSL193"  "pM11003" "pM4175"
# good Reference_Pair is the variable for CNV
# we need to map Participant_ID to Reference_Pair
# mm.meta$`Sample ID`[mm.meta$`Sample ID`%nin% TP53$ID] %in% clinical.and.terra.ref.pairs.samples$Participant_ID

TP53 <- TP53 |> 
  left_join(clinical.and.terra.ref.pairs.samples |> select(Participant_ID, Reference_Pair), by=c("ID"="Reference_Pair")) |>
  mutate(`Sample ID` = case_when(dataset == "SU2C" ~  Participant_ID,
                                 TRUE ~ ID))

TP53 <- TP53 %>% mutate(across(where(is.numeric), ~{1*(.x>.14) + 1*(.x>1.73)}))

mm <- mm %>% left_join(TP53 %>% select(`Sample ID`, TP53_CN), by="Sample ID") |> rename(ID=`Sample ID`)

#mm is missing TP53 status for NDMM from SU2C
mm |> filter(is.na(TP53_CN))

write_tsv(x = data.frame(mm), file = "../data/tmp_rescore.tsv")


# Score -------------------------------------------------------------------

mm <- mm |> 
  mutate(Sex=Gender, Disease_Status=Stage, cohort=Cohort, Days=NA, Progression_status=NA)

sig.su2c.jb <- function(df){
  df$t14_16 <- ifelse(df$t14_16 >=1, -1, 0)
  df$t14_20 <- ifelse(df$t14_20 >=1, -1, 0) 
  df$KRAS <- ifelse(df$SNV_KRAS >=1, 1, 0)
  df$NRAS <- ifelse(df$SNV_NRAS >=1, 1, 0)
  df$FAM46C <- ifelse(df$SNV_FAM46C >=1, 1, 0)
  df$HyperDiploidy <- ifelse(df$HyperDiploidy ==1, 1, 0)
  df$del_chr_1p <- ifelse(df$del_G.peak.6_del_1p22.1 >=1 | df$del_G.peak.7_del_1p12 >=1 | df$del_G.peak.5_del_1p32.3 >=1 , 1, 0)
  df$amp_chr_1q <- ifelse(df$amp_chr_1q >=1 | df$amp_G.peak.1_amp_1q21.2 >=1, 1, 0)
  df$amp_chr_3q <- ifelse(df$amp_G.peak.2_amp_3q26.2 >=1 & df$HyperDiploidy==0, 1, 0) # don't double count HRD
  df$del_chr_4p <- ifelse(df$del_G.peak.11_del_4p16.3 >=1, 1, 0)
  df$del_chr_4q <- ifelse(df$del_G.peak.12_del_4q34.3 >=1, 1, 0)
  df$del_chr_8p <- ifelse(df$del_chr_8p >=1 | df$del_G.peak.16_del_8p23.3 >=1, 1, 0)
  df$amp_chr_8q <- ifelse(df$del_G.peak.17_del_8q24.21 >=1 | df$del_G.peak.18_del_8q24.21 >=1 | df$amp_G.peak.3_amp_8q24.21 >=1 | df$MYC_translocation >=1, 1, 0)
  df$del_chr_12p <- ifelse(df$del_G.peak.22_del_12p13.33 >=1, 1, 0)
  df$del_chr_14q <- ifelse(df$del_G.peak.27_del_14q23.3 >=1 | df$del_G.peak.28_del_14q31.1 >=1, 1, 0)
  df$del_chr_16q <- ifelse(df$del_G.peak.30_del_16q12.1 >=1 | df$del_chr_16q >=1, 1, 0)
  df$del_chr_17q <- ifelse(df$del_G.peak.31_del_17q21.2 >=1, 1, 0)
  df$total_sig <- rowSums(df[,c("t14_16", "t14_20", "KRAS", "NRAS", "FAM46C", "HyperDiploidy", "del_chr_1p", "amp_chr_1q", "amp_chr_3q", "del_chr_4p", "del_chr_4q", "del_chr_8p", "amp_chr_8q", "del_chr_12p", "del_chr_14q", "del_chr_16q", "del_chr_17q")])
  df %>% select(ID, Sex, Age, Disease_Status, cohort, t14_16, t14_20, KRAS, NRAS, FAM46C, HyperDiploidy, del_chr_1p, amp_chr_1q, amp_chr_3q, del_chr_4p, del_chr_4q, del_chr_8p, amp_chr_8q, del_chr_12p, del_chr_14q, del_chr_16q, del_chr_17q, total_sig, Days, Progression_status)
}


# Fig 1B -----------
sig.mm <- sig.su2c.jb(mm)

sig.mm |> write_tsv("../data/export_MM-like_score.tsv")

sig.mm |> dunn_test(total_sig ~ Disease_Status)
sig.mm |> kruskal_test(total_sig ~ Disease_Status)

sig.mm.stats <- sig.mm |> group_by(Disease_Status) |> summarise(n=n(), label=paste0("N=", n), y=12)

fig1b <- ggplot(sig.mm, aes(factor(Disease_Status, levels=c("MGUS", "SMM", "MM")), total_sig, fill=factor(Disease_Status, levels=c("MGUS", "SMM", "MM")))) + 
  geom_violin(scale = "width", adjust=1.5, alpha=0.25) +
  geom_boxplot(width=0.4) +
  geom_text(data=sig.mm.stats, aes(y=y, label=label), size=2.5) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option = "G") +
  scale_y_continuous(breaks=seq(-1, 11)) +
  labs(x="", y="MM-like score", fill="") + 
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none")

print(fig1b)

ggsave(filename = "../Figures/Signature_Disease_stages.pdf", plot = fig1b, width = 2, height = 2)
ggsave(filename = "../Figures/Signature_Disease_stages.png", plot = fig1b, width = 2, height = 2)

# 
# SU2C.sig <- read_csv("../JB/Modified dataframes/Full_Signature_SU2C_only.csv")
# SU2C.Mark.sig <- read_csv("../JB/Modified dataframes/SU2C_Mark_Signature_SU2C_only.csv")
SU2C.Mark.survival <- readxl::read_xlsx("../JB/Original data/Mark.progression.targets.87.updated.xlsx")

# SU2C.sig |> 
#   filter(Disease_Status%in%c("Low", "Intermediate", "High")) |>
#   View()

# Fig 1C ---------
sig.smm <- sig.mm |> filter(Disease_Status=="SMM") |> left_join(mm |> select(1:9))
table(sig.smm$cohort, sig.smm$StageAnd22020)
# problem: 61 bustoros don't have stage
# stage is in SU2C_Mark_Signature_SU2C_only.csv

sig.smm.non.null <- sig.smm |> 
  left_join(SU2C.Mark.survival |> select(`Sample ID`, RISK), by=c("ID"="Sample ID")) |>
  mutate(StageAnd22020 = case_when(StageAnd22020 == "SMM" & !is.na(RISK) ~ RISK,
                                   TRUE ~ StageAnd22020)) |>
  mutate(risk_class=case_when(
    StageAnd22020 %in% c("HRSMM", "High") ~ "High",
    StageAnd22020 %in% c("LRSMM", "Low") ~ "Low",
    StageAnd22020 %in% c("IRSMM", "Intermediate") ~ "Intermediate",
    TRUE ~ "SMM")) |>
  filter(risk_class!="SMM")
  
# sig.smm$ID[sig.smm$StageAnd22020=="SMM"] %in% SU2C.sig$ID
# rescue Mark

sig.smm.non.null |> write_tsv("../data/export_MM-like_score_SMM.tsv")

sig.smm.non.null |> dunn_test(total_sig ~ risk_class)
sig.smm.non.null |> kruskal_test(total_sig ~ risk_class)

sig.smm.stats <- sig.smm.non.null |> group_by(risk_class) |> summarise(n=n(), label=paste0("N=", n), y=9)

fig1c <- ggplot(sig.smm.non.null, aes(factor(risk_class, levels=c("Low", "Intermediate", "High")), total_sig, fill=factor(risk_class, levels=c("Low", "Intermediate", "High")))) + 
  geom_violin(scale = "width", adjust=1.5, alpha=0.25) +
  geom_boxplot(width=0.4) +
  geom_text(data=sig.smm.stats, aes(y=y, label=label), size=2.5) +
  scale_fill_brewer(palette = 1) +
  scale_y_continuous(breaks=seq(-1, 11)) +
  labs(x="20/2/20 Risk", y="MM-like score", fill="") + 
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none")

print(fig1c)

ggsave(filename = "../figures/SMM_Risk.pdf", plot=fig1c, width = 2, height = 2)
ggsave(filename = "../figures/SMM_Risk.png", plot=fig1c, width = 2, height = 2)


# Fig 1D ----------

# where can I find survival data?
survival.su2c <- readxl::read_xlsx("../JB/Original data/SU2C.clinical.data.xlsx")

# summary(coxph(Surv(Days, Progression_status) ~ total_sig, data=sig.su2c.s))

# > survival.su2c$Participant_ID_JB[survival.su2c$Participant_ID_JB%nin% sig.mm$ID] 
# [1] "pM4175"  "pM11003"
# great

survival.sig <- sig.smm |> 
  filter(Disease_Status=="SMM") |>
  left_join(SU2C.Mark.survival |> select(`Sample ID`, RISK), by=c("ID"="Sample ID")) |>
  mutate(StageAnd22020 = case_when(StageAnd22020 == "SMM" & !is.na(RISK) ~ RISK,
                                   TRUE ~ StageAnd22020)) |>
  mutate(risk_class=case_when(
    StageAnd22020 %in% c("HRSMM", "High") ~ "High",
    StageAnd22020 %in% c("LRSMM", "Low") ~ "Low",
    StageAnd22020 %in% c("IRSMM", "Intermediate") ~ "Intermediate",
    TRUE ~ "SMM")) |>
  select(-Days, -Progression_status) |>
  inner_join(survival.su2c |> select(Participant_ID_JB, `Days (PFS sample date, treatment censor)`, Progression), by=c("ID"="Participant_ID_JB")) |>
  rename(Days=`Days (PFS sample date, treatment censor)`) |>
  rename(Progression_status=Progression) |>
  filter(Days>=90)
  
survival.sig |> write_tsv("../data/export_MM-like_score_SMM_survival_SU2C.tsv")

# survival.su2c  |> filter(`Days (PFS sample date, treatment censor)`<90) |> View()
dim(survival.sig)

summary(coxph(Surv(Days, Progression_status) ~ total_sig, data=survival.sig))
summary(coxph(Surv(Days, Progression_status) ~ total_sig>1, data=survival.sig))
summary(coxph(Surv(Days, Progression_status) ~ total_sig+StageAnd22020, data=survival.sig))
# add staging

survival.sig$Score = survival.sig$total_sig >= 2

fit <- survfit(Surv(Days, Progression_status) ~ Score, data=survival.sig)
summary(coxph(Surv(Days, Progression_status) ~ Score, data=survival.sig))
exp(confint(coxph(Surv(Days, Progression_status) ~ Score, data=survival.sig)))

surv.plot <- ggsurvplot(fit, 
           data = survival.sig,
           surv.median.line = "none", # Add medians survival
           
           # Change legends: title & labels
           # Add p-value and interval
           pval = FALSE,
           conf.int = TRUE,
           xscale = "d_y",
           ylab = c("Probability of Progression"),
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.25,
           tables.theme = theme_cleantable(),
           xlim = c(0, 365.25*5),
           break.time.by = 365.25,
           fun = "event",
           
           # strata group labels
           legend.labs = c("Score <= 1", "Score > 1"),
           
           # Color palettes
           palette = c("#377EB8", "#E41A1C"),
           ggtheme = theme_bw(7), # Change ggplot2 theme
           fontsize = 2.5)

#thanks to https://github.com/kassambara/survminer/issues/152#issuecomment-320886179
pdf("../figures/km.s.plot.pdf", width = 2.4, height = 2.4)
print(surv.plot, newpage = FALSE)
dev.off()
