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
  df %>% select(`Sample ID`, Sex, Age, Disease_Status, cohort, t14_16, t14_20, KRAS, NRAS, FAM46C, HyperDiploidy, del_chr_1p, amp_chr_1q, amp_chr_3q, del_chr_4p, del_chr_4q, del_chr_8p, amp_chr_8q, del_chr_12p, del_chr_14q, del_chr_16q, del_chr_17q, total_sig, Days, Progression_status)
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
  geom_text(data=sig.mm.stats, aes(y=y, label=label), size=6/.pt) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option = "G") +
  scale_y_continuous(breaks=seq(-1, 11)) +
  labs(x="", y="MM-like score", fill="") + 
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text = element_text(size = 6))

print(fig1b)

ggsave(filename = "../Figures/Signature_Disease_stages.pdf", plot = fig1b, width = 2, height = 2)
ggsave(filename = "../Figures/Signature_Disease_stages.png", plot = fig1b, width = 2, height = 2)

SU2C.Mark.survival <- readxl::read_xlsx("../JB/Original data/Mark.progression.targets.87.updated.xlsx")

# Fig 1C ---------
sig.smm <- sig.mm |> filter(Disease_Status=="SMM") |> left_join(mm |> select(1:9))
table(sig.smm$cohort, sig.smm$StageAnd22020)
# problem: 61 bustoros don't have stage
# stage is in SU2C_Mark_Signature_SU2C_only.csv

sig.smm.non.null <- sig.smm |> 
  left_join(SU2C.Mark.survival |> select(`Sample ID`, RISK), by=c("Sample ID"="Sample ID")) |>
  mutate(StageAnd22020 = case_when(StageAnd22020 == "SMM" & !is.na(RISK) ~ RISK,
                                   TRUE ~ StageAnd22020)) |>
  mutate(risk_class=case_when(
    StageAnd22020 %in% c("HRSMM", "High") ~ "High",
    StageAnd22020 %in% c("LRSMM", "Low") ~ "Low",
    StageAnd22020 %in% c("IRSMM", "Intermediate") ~ "Intermediate",
    TRUE ~ "SMM")) |>
  filter(risk_class!="SMM")

sig.smm.non.null |> write_tsv("../data/export_MM-like_score_SMM.tsv")

sig.smm.non.null |> dunn_test(total_sig ~ risk_class)
# # A tibble: 3 × 9
# .y.       group1       group2          n1    n2 statistic          p     p.adj p.adj.signif
# * <chr>     <chr>        <chr>        <int> <int>     <dbl>      <dbl>     <dbl> <chr>       
#   1 total_sig High         Intermediate    64    32     -2.41 0.0161     0.0321    *           
#   2 total_sig High         Low             64    54     -4.61 0.00000412 0.0000124 ****        
#   3 total_sig Intermediate Low             32    54     -1.48 0.139      0.139     ns    
sig.smm.non.null |> kruskal_test(total_sig ~ risk_class)
# # A tibble: 1 × 6
# .y.           n statistic    df         p method        
# * <chr>     <int>     <dbl> <int>     <dbl> <chr>         
#   1 total_sig   150      21.6     2 0.0000199 Kruskal-Wallis

sig.smm.stats <- sig.smm.non.null |> group_by(risk_class) |> summarise(n=n(), label=paste0("N=", n), y=9)

fig1c <- ggplot(sig.smm.non.null, aes(factor(risk_class, levels=c("Low", "Intermediate", "High")), total_sig, fill=factor(risk_class, levels=c("Low", "Intermediate", "High")))) + 
  geom_violin(scale = "width", adjust=1.5, alpha=0.25) +
  geom_boxplot(width=0.4) +
  geom_text(data=sig.smm.stats, aes(y=y, label=label), size=6/.pt) +
  scale_fill_brewer(palette = 1) +
  scale_y_continuous(breaks=seq(-1, 11)) +
  labs(x="20/2/20 Risk", y="MM-like score", fill="") + 
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text = element_text(size=6))

print(fig1c)

ggsave(filename = "../figures/SMM_Risk.pdf", plot=fig1c, width = 2, height = 2)
ggsave(filename = "../figures/SMM_Risk.png", plot=fig1c, width = 2, height = 2)


# Fig 1D ----------

survival.su2c <- readxl::read_xlsx("../JB/Original data/SU2C.clinical.data.xlsx")

survival.sig <- sig.smm |> 
  filter(Disease_Status=="SMM") |>
  left_join(SU2C.Mark.survival |> select(`Sample ID`, RISK), by=c("Sample ID"="Sample ID")) |>
  mutate(StageAnd22020 = case_when(StageAnd22020 == "SMM" & !is.na(RISK) ~ RISK,
                                   TRUE ~ StageAnd22020)) |>
  mutate(risk_class=case_when(
    StageAnd22020 %in% c("HRSMM", "High") ~ "High",
    StageAnd22020 %in% c("LRSMM", "Low") ~ "Low",
    StageAnd22020 %in% c("IRSMM", "Intermediate") ~ "Intermediate",
    TRUE ~ "SMM")) |>
  select(-Days, -Progression_status) |>
  inner_join(survival.su2c |> select(Participant_ID_JB, `Days (PFS sample date, treatment censor)`, Progression), by=c("Sample ID"="Participant_ID_JB")) |>
  rename(Days=`Days (PFS sample date, treatment censor)`) |>
  rename(Progression_status=Progression) |>
  filter(Days>=90)

survival.sig |> write_tsv("../data/export_MM-like_score_SMM_survival_SU2C.tsv")

# survival.su2c  |> filter(`Days (PFS sample date, treatment censor)`<90) |> View()
dim(survival.sig)

survival.sig$HRSMM=survival.sig$StageAnd22020=="HRSMM"

summary(coxph(Surv(Days, Progression_status) ~ total_sig, data=survival.sig))
# to limit confusion, we will report the dichotomized p- value together with HR. 
# This is to stay consistent with the KM plot which has two grups.
# In the future this should probably be interpreted as a "per point" signature.
summary(coxph(Surv(Days, Progression_status) ~ total_sig>1, data=survival.sig))
# Call:
#   coxph(formula = Surv(Days, Progression_status) ~ total_sig > 
#           1, data = survival.sig)
# 
# n= 63, number of events= 13 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# total_sig > 1TRUE 1.4514    4.2692   0.6681 2.172   0.0298 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# total_sig > 1TRUE     4.269     0.2342     1.153     15.81
# 
# Concordance= 0.622  (se = 0.086 )
# Likelihood ratio test= 5.6  on 1 df,   p=0.02
# Wald test            = 4.72  on 1 df,   p=0.03
# Score (logrank) test = 5.59  on 1 df,   p=0.02


# Fig 1E (from Sofia) -----------------------------------------------------

# KM - Signature Bustoros only
m.sig.su2c <- read_csv("../JB/Modified dataframes/Mark_Signature_SU2C_only.csv")
m.sig.su2c <-  m.sig.su2c %>% mutate(Disease_Status=factor(Disease_Status, levels = c("Low", "Intermediate", "High")))
m.sig.su2c <-  m.sig.su2c %>% mutate(HRSMM=Disease_Status=="High")

modify.score <- function(data) {
  
  #dichotomize
  data <- data %>% 
    mutate(Score = as.numeric(total_sig>=2))
}

km.m <- modify.score(m.sig.su2c)

table(km.m$Score)

km.m$Score=as.logical(km.m$Score)
summary(coxph(Surv(Days, Progression_status) ~ Score, data=km.m))

# Call:
#   coxph(formula = Surv(Days, Progression_status) ~ Score, data = km.m)
# 
# n= 87, number of events= 58 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# ScoreTRUE 1.2218    3.3935   0.2807 4.353 1.34e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# ScoreTRUE     3.393     0.2947     1.958     5.883
# 
# Concordance= 0.653  (se = 0.031 )
# Likelihood ratio test= 19.5  on 1 df,   p=1e-05
# Wald test            = 18.95  on 1 df,   p=1e-05
# Score (logrank) test = 21.07  on 1 df,   p=4e-06

fit_Score_Bustoros <- coxph(Surv(Days, Progression_status) ~ Score+Disease_Status, data=as.data.frame(km.m))
fit_Score_Bustoros

ggforest.jco <- ggforest(fit_Score_Bustoros, fontsize = 1, noDigits = 2,main = NULL)

ggsave(filename = "../figures/forest-jco.pdf", ggforest.jco, width = 2.5, height = 2)


#KM
km <- function(data){
  
  fit <- survfit(Surv(Days, Progression_status) ~ Score, data=data)
  
  ggsurvplot(fit, 
             data = data,
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
             ggtheme = theme_bw(16), # Change ggplot2 theme
             fontsize = 5)
}

km(km.m)

# Revisions ---------------------------------------------------------------

# reviewer 1, is it significant in multivariate
survival.sig$StageAnd22020
survival.sig$StageAnd22020 <- factor(survival.sig$StageAnd22020, levels=c("LRSMM", "IRSMM", "HRSMM"))
summary(coxph(Surv(Days, Progression_status) ~ HRSMM, data=survival.sig))
summary(coxph(Surv(Days, Progression_status) ~ StageAnd22020, data=survival.sig))
# add staging

survival.sig$Score = survival.sig$total_sig >= 2

fit <- survfit(Surv(Days, Progression_status) ~ Score, data=survival.sig)
summary(coxph(Surv(Days, Progression_status) ~ Score, data=survival.sig))
summary(coxph(Surv(Days, Progression_status) ~ Score+HRSMM, data=survival.sig))
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

pdf("../figures/km.s.plot.pdf", width = 2.4, height = 2.4)
print(surv.plot, newpage = FALSE)
dev.off()

# Extended Data Figure 2D -------------------------------------------------

# thanks to https://stackoverflow.com/questions/75152731/r-ggforest-error-in-match-namesclabs-namesxi-names-do-not-match-previou

fit_Score_SU2C <- coxph(Surv(Days, Progression_status) ~ Score+StageAnd22020, data=as.data.frame(survival.sig))
summary(fit_Score_SU2C)
ggforest.su2c <- ggforest(fit_Score_SU2C, fontsize = 1, noDigits = 2,main = NULL, refLabel = c("LR"))

ggsave(filename = "../figures/forest-su2c.pdf", , width = 2.5, height = 2)
