setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(RColorBrewer)

# Merge genomics and survival ---------------------------------------------

xlsxfile <- "../data/Boyle_NatComm_SMM_data.xlsx"
excel_sheets(xlsxfile)

smm89 <- read_xlsx(xlsxfile, sheet = "SMM89")
colnames(smm89)
smm89survival <- read_xlsx(xlsxfile, sheet = "Survival")

table(smm89$pat_id %in% smm89survival$fnpatid)

smm89m <- smm89 |> inner_join(smm89survival, by=c("pat_id"="fnpatid"))

sig.smm89.jb <- function(df){
  df$t14_16 <- ifelse(df$transloc_maf >=1, -1, 0)
  df$t14_20 <- ifelse(df$transloc_mafb >=1, -1, 0) 
  df$KRAS <- ifelse(df$mutation_kras >=1, 1, 0)
  df$NRAS <- ifelse(df$mutation_nras >=1, 1, 0)
  df$FAM46C <- ifelse(df$mutation_fam46c >=1, 1, 0)
  df$HyperDiploidy <- ifelse(df$hrd ==1, 1, 0)
  df$del_chr_1p <- ifelse(df$copynum_cdkn2c <= 1 | df$copynum_fam46c <= 1 , 1, 0)
  df$amp_chr_1q <- ifelse(df$copynum_1q21_3 >=3, 1, 0)
  # df$amp_chr_3q <- ifelse(df$amp_G.peak.2_amp_3q26.2 >=1 & df$HyperDiploidy==0, 1, 0) # don't double count HRD
  df$del_chr_4p <- ifelse(df$copynum_whsc1l <=1, 1, 0) # corresponds to 4p16.3 indeed
  # df$del_chr_4q <- ifelse(df$del_G.peak.12_del_4q34.3 >=1, 1, 0) # no gene...
  # df$del_chr_8p <- ifelse(df$del_chr_8p >=1 | df$del_G.peak.16_del_8p23.3 >=1, 1, 0) # no gene ...
  df$amp_chr_8q <- ifelse(df$copynum_myc >=3 | df$transloc_myc >=1, 1, 0)
  df$del_chr_12p <- ifelse(df$copynum_cdkn1b <=1, 1, 0)
  df$del_chr_14q <- ifelse(df$copynum_traf3 <=1, 1, 0)
  df$del_chr_16q <- ifelse(df$copynum_cyld <=1, 1, 0)
  # df$del_chr_17q <- ifelse(df$del_G.peak.31_del_17q21.2 >=1, 1, 0)
  # df$total_sig <- rowSums(df[,c("t14_16", "t14_20", "KRAS", "NRAS", "FAM46C", "HyperDiploidy", "del_chr_1p", "amp_chr_1q", "amp_chr_3q", "del_chr_4p", "del_chr_4q", "del_chr_8p", "amp_chr_8q", "del_chr_12p", "del_chr_14q", "del_chr_16q", "del_chr_17q")])
  df$total_sig <- rowSums(df[,c("t14_16", "t14_20", "KRAS", "NRAS", "FAM46C", "HyperDiploidy", "del_chr_1p", "amp_chr_1q", "del_chr_4p", "amp_chr_8q", "del_chr_12p", "del_chr_14q", "del_chr_16q")])
  df
}

sig.smm89.jb.noCN <- function(df){
  df$t14_16 <- ifelse(df$transloc_maf >=1, -1, 0)
  df$t14_20 <- ifelse(df$transloc_mafb >=1, -1, 0) 
  df$KRAS <- ifelse(df$mutation_kras >=1, 1, 0)
  df$NRAS <- ifelse(df$mutation_nras >=1, 1, 0)
  df$FAM46C <- ifelse(df$mutation_fam46c >=1, 1, 0)
  df$HyperDiploidy <- ifelse(df$hrd ==1, 1, 0)
  df$total_sig_noCN <- rowSums(df[,c("t14_16", "t14_20", "KRAS", "NRAS", "FAM46C", "HyperDiploidy")])
  df
}

smm89m <- sig.smm89.jb(smm89m)
smm89m <- sig.smm89.jb.noCN(smm89m)

smm89m$sig_class_noCN <- ifelse(smm89m$total_sig_noCN >= 2, ">1", "≤1")
smm89m$sig_class <- ifelse(smm89m$total_sig >= 2, ">1", "≤1")
smm89m$S = smm89m$total_sig_noCN
smm89m$Score = smm89m$total_sig_noCN>1


smm89m |> 
  select(main_id,pat_id,Mayo,PFS_time,PFS_status,total_sig,sig_class,total_sig_noCN,sig_class_noCN,S,Score) |> 
  write_tsv("../data/Export_Boyle_with_score.tsv")

## with High Risk
smm89m.risk <- smm89m |> 
  filter(Mayo!="NA") |> 
  mutate(HRSMM=Mayo=="2")
smm89m.risk$Mayo = factor(smm89m.risk$Mayo, levels=c("0", "1", "2"))

smm89m.risk |> 
  select(main_id,pat_id,Mayo,PFS_time,PFS_status,total_sig,sig_class,total_sig_noCN,sig_class_noCN,S,Score) |> 
  write_tsv("../data/Export_Boyle_Mayo_risk_with_score.tsv")


# Cox model ---------------------------------------------------------------

smm89m$sig_class_noCN <- factor(smm89m$sig_class_noCN, levels=c("≤1", ">1"))
smm89m$sig_class <- factor(smm89m$sig_class, levels=c("≤1", ">1"))

## CN ----
summary(coxph(Surv(PFS_time, PFS_status) ~ total_sig, data=smm89m))
exp(confint(coxph(Surv(PFS_time, PFS_status) ~ total_sig, data=smm89m)))
# HR: 1.4137, CI 95% 1.023987 1.951607, P=0.03

## No CN ----
summary(coxph(Surv(PFS_time, PFS_status) ~ total_sig_noCN, data=smm89m))
exp(confint(coxph(Surv(PFS_time, PFS_status) ~ total_sig_noCN, data=smm89m)))
# HR: 1.7732, CI 95% 1.048593 2.998673, P=0.03

summary(coxph(Surv(PFS_time, PFS_status) ~ total_sig + Mayo, data=smm89m |> filter(Mayo!="NA") |> mutate(HRSMM=Mayo=="2")))
summary(coxph(Surv(PFS_time, PFS_status) ~ total_sig_noCN + HRSMM, data=smm89m |> filter(Mayo!="NA") |> mutate(HRSMM=Mayo=="2")))
summary(coxph(Surv(PFS_time, PFS_status) ~ sig_class_noCN + HRSMM, data=smm89m |> filter(Mayo!="NA") |> mutate(HRSMM=Mayo=="2")))

# 20/2/20 Risk ----
summary(coxph(Surv(PFS_time, PFS_status) ~ Mayo, data=smm89m.risk))
summary(coxph(Surv(PFS_time, PFS_status) ~ HRSMM, data=smm89m.risk))

# Figure ------------------------------------------------------------------

## no CN ----

summary(coxph(Surv(PFS_time, PFS_status) ~ total_sig_noCN, data=smm89m|> filter(PFS_time>0)))
exp(confint(coxph(Surv(PFS_time, PFS_status) ~ total_sig_noCN, data=smm89m|> filter(PFS_time>0))))
fit_noCN <- survfit(Surv(PFS_time, PFS_status) ~ S, data=smm89m |> filter(PFS_time>0))

surv.plot <- ggsurvplot(fit_noCN,
                        data = smm89m,
                        surv.median.line = "none", # Add medians survival
                        
                        # Change legends: title & labels
                        # Add p-value and interval
                        pval = FALSE,
                        conf.int = FALSE,
                        # xscale = "d_y",
                        ylab = c("Probability of Progression"),
                        # Add risk table
                        risk.table = TRUE,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        # xlim = c(0, 365.25*5),
                        # break.time.by = 365.25,
                        fun = "event",
                        
                        # strata group labels
                        legend.labs = as.character(seq(-1,2)),
                        
                        # Color palettes
                        palette = brewer.pal(n = 6, name = "RdBu")[c(6, 5, 4, 2)],
                        ggtheme = theme_bw(7), # Change ggplot2 theme
                        fontsize = 2)

pdf("../figures/km.smm89.boyle.plot.all.scores.pdf", width = 2.4, height = 2.4)
print(surv.plot, newpage = FALSE)
dev.off()

# Forest plot -------------------------------------------------------------

fit_noCN_HRSMM <- coxph(Surv(PFS_time, PFS_status) ~ Score + Mayo, data=as.data.frame(smm89m.risk))
# summary(fit_noCN_HRSMM)
# Call:
#   coxph(formula = Surv(PFS_time, PFS_status) ~ Score + Mayo, data = as.data.frame(smm89m.risk))
# 
# n= 75, number of events= 26 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# ScoreTRUE 0.9538    2.5956   0.4765 2.002 0.045319 *  
#   Mayo1     0.5061    1.6589   0.5476 0.924 0.355361    
# Mayo2     1.6858    5.3966   0.4909 3.434 0.000595 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# ScoreTRUE     2.596     0.3853    1.0201     6.604
# Mayo1         1.659     0.6028    0.5671     4.852
# Mayo2         5.397     0.1853    2.0619    14.124
# 
# Concordance= 0.763  (se = 0.042 )
# Likelihood ratio test= 15.54  on 3 df,   p=0.001
# Wald test            = 16.81  on 3 df,   p=8e-04
# Score (logrank) test = 19.86  on 3 df,   p=2e-04
ggforest.boyle <- ggforest(fit_noCN_HRSMM, fontsize = 1, noDigits = 2,main = NULL)

ggsave(filename = "../figures/forest-boyle.pdf", , width = 2.5, height = 2)

