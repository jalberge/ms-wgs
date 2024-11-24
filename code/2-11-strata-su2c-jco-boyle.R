setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(RColorBrewer)
library(broom)

# Strata based MV ---------------------------------------------------------

## Boyle et al. ----

smm89m <- read_tsv("../data/Export_Boyle_with_score.tsv")
smm89m$S = smm89m$total_sig_noCN
smm89m$Score = smm89m$total_sig_noCN>1
smm89m.risk <- smm89m |> filter(Mayo!="NA")
smm89m.risk$IMWG_risk <- factor(smm89m.risk$Mayo, levels=c(0, 1, 2), labels=c("Low", "Intermediate", "High"))
smm89m.risk$days <- smm89m.risk$PFS_time*365.25
smm89m.risk$event <- smm89m.risk$PFS_status
smm89m.risk$study <- "NatComm"

## SU2C / JCO ----

# load score from SMM
smm.scores <- read_tsv("../data/20241110_MM-like_IMWG_cytogenet_scores.tsv")

# survival SU2C
survival.SU2C <- read_tsv("../data/export_MM-like_score_SMM_survival_SU2C.tsv")
survival.SU2C$ID=survival.SU2C$`Sample ID`
# survival MARK
survival.mark <- read_csv("../JB/Modified dataframes/Mark_Signature_SU2C_only.csv")

mark.su2c.survival <- full_join(survival.SU2C, survival.mark)
mark.su2c.survival$study <- ifelse(mark.su2c.survival$cohort%in% c("SU2C", "CTC"), "SU2C", "JCO")


# Combine ---------------------------------------------------------

mark.su2c.survival$days <- mark.su2c.survival$Days
mark.su2c.survival$event <- mark.su2c.survival$Progression_status
mark.su2c.survival$Score <- mark.su2c.survival$total_sig > 1

mark.su2c.survival <- mark.su2c.survival |> 
  mutate(IMWG_risk=case_when(study=="SU2C" ~ risk_class, study=="JCO" ~ Disease_Status))

smm89m.risk.simple <- smm89m.risk |> select(days, event, Score, IMWG_risk, study)
mark.su2c.survival.simple <- mark.su2c.survival |> select(days, event, Score, IMWG_risk, study)
mark.su2c.boyle <- bind_rows(smm89m.risk.simple, mark.su2c.survival.simple)

# Strata analysis -------------------------------------------------

# mark.su2c.boyle$study <- as.numeric(factor(mark.su2c.boyle$study))
mark.su2c.boyle$IMWG_risk <- factor(mark.su2c.boyle$IMWG_risk, levels=c("Low", "Intermediate", "High"))
strata.model <- coxph(Surv(days, event) ~ Score + IMWG_risk + strata(study), data=as.data.frame(mark.su2c.boyle))
summary(strata.model)

# n= 225, number of events= 97 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# ScoreTRUE             1.0311    2.8040   0.2304 4.476 7.61e-06 ***
#   IMWG_riskIntermediate 0.7658    2.1508   0.3352 2.285   0.0223 *  
#   IMWG_riskHigh         1.8034    6.0700   0.3021 5.969 2.39e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# ScoreTRUE                 2.804     0.3566     1.785     4.404
# IMWG_riskIntermediate     2.151     0.4650     1.115     4.149
# IMWG_riskHigh             6.070     0.1647     3.357    10.974
# 
# Concordance= 0.77  (se = 0.024 )
# Likelihood ratio test= 73.81  on 3 df,   p=7e-16
# Wald test            = 66.55  on 3 df,   p=2e-14
# Score (logrank) test = 79.33  on 3 df,   p=<2e-16


strata.df <- tidy(strata.model, exponentiate = TRUE, conf.int = TRUE)
strata.df[strata.df["term"]=="ScoreTRUE", "term"] <- "Score > 1"
strata.df[strata.df["term"]=="IMWG_riskIntermediate", "term"] <- "Intermediate"
strata.df[strata.df["term"]=="IMWG_riskHigh", "term"] <- "High"

# 3 patients have no follow up (removed from the model above - here for stats)
mark.su2c.boyle |> filter(days>1) |> group_by(IMWG_risk) |> count()
mark.su2c.boyle |> filter(days>1) |> group_by(study) |> count()
mark.su2c.boyle |> filter(days>1) |> group_by(study, IMWG_risk) |> count()

stats.df <- as.data.frame(glance(strata.model))

stats.text <- paste0("# Events: ", stats.df["nevent"], "; Global p-value (Log-Rank): ", format(stats.df["p.value.log"], digits=2), 
       "\n", 
       "AIC: ", format(stats.df["AIC"], digits=4), "; Concordance Index: ", format(stats.df["concordance"], digits=1))

stratified.survival <- ggplot(strata.df, aes(y=term)) +
  geom_vline(xintercept = 1, linetype=2, color="grey") +
  geom_point(aes(x=estimate), shape = 15, size=4) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0.1) +
  geom_text(aes(x=8, label=format(p.value, digits=2)), position = position_nudge(y=0.2), fontface = "italic") +
  labs(x=stats.text, y="") +
  scale_x_continuous(transform = "pseudo_log", limits=c(0, NA)) +
  theme_classic2()

ggsave("../figures/stratified.survival.pdf", plot = stratified.survival, width = 3, height = 3, scale = 2)
ggsave("../figures/stratified.survival.png", plot = stratified.survival, width = 3, height = 3, scale = 2)
# bug with strata
# 
# ggforest.strata <- ggforest(strata.model, fontsize = 1, noDigits = 2,main = NULL)
# 
# ggsave(filename = "../figures/forest-strata.pdf", , width = 2.5, height = 2)
