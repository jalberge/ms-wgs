setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(survival)
library(viridis)
library(survminer)

source("0_annotate_samples.R")


# MM-like score from SMM --------------------------------------------------

# load score from SMM
sig.smm.non.null <- read_tsv("../data/export_MM-like_score_SMM.tsv")

# IMWG risk ----
# add IMWG High Risk
muts.mm <- read_tsv("../data/fig1_all_matrix_patients_with_disease_stage.tsv")
muts.mm$IMWG_FISH_HighRisk <- 1 * as.numeric( (muts.mm$`t(14;16)(MAF)` + muts.mm$`t(4;14)(MMSET)` + 
                                                 muts.mm$`-13q` + muts.mm$`-13q14.3` + muts.mm$`-13q31.2` + muts.mm$`-13q12.11` +
                                                 muts.mm$`+1q`+muts.mm$`+1q21.2`+muts.mm$`+1q21.2`)>=1 )
all( sig.smm.non.null$ID %in% muts.mm$Participant_ID )

sig.smm.non.null <- sig.smm.non.null |> left_join(muts.mm |> select(Participant_ID, IMWG_FISH_HighRisk), by=c("ID"="Participant_ID"))

# in order to re-assess score, we need exact number of 20/2/2points


## From Bustoros -----------------------------------------------------------

# first; jco paper data
all( clinical.annotation$Participant_ID %in% muts.mm$Participant_ID )
sig.smm.non.null
# BM, FLC and M spike for Mark are here
SU2C.Mark.survival <- readxl::read_xlsx("../JB/Original data/Mark.progression.targets.87.updated.xlsx")
SU2C.Scores.For.FISH.IMWG <- SU2C.Mark.survival |> 
  # filter(RISK=="High") |> 
  select(`Sample ID`, RISK, BM_Involvement_per_Bx, M_spike, Inv_Uninv_LC_Ratio, t14_16, del13q, gain1q, t4_14) |> 
  rowwise() |>
  mutate(HRFISH=
           as.numeric(sum(as.numeric(t14_16), as.numeric(del13q), as.numeric(gain1q), as.numeric(t4_14), na.rm=TRUE)>=1)) |> 
  mutate(`20_2_20_points`=sum(BM_Involvement_per_Bx>0.20, M_spike>2, Inv_Uninv_LC_Ratio>20)) |> 
  mutate(IMWG_FISH_points=`20_2_20_points`+HRFISH) |>
  mutate(IMWG_FISH_class=case_when(
    IMWG_FISH_points==0 ~ "Low",
    IMWG_FISH_points==1 ~ "LowInt",
    IMWG_FISH_points==2 ~ "Int",
    IMWG_FISH_points>=3 ~ "High",
    TRUE ~ "NA"
  ))

all( pull( sig.smm.non.null[sig.smm.non.null$cohort=="Bustoros",], "ID") %in% SU2C.Scores.For.FISH.IMWG$`Sample ID` )
# TRUE
sig.smm.non.null <- sig.smm.non.null |> left_join(SU2C.Scores.For.FISH.IMWG |> select(`Sample ID`, `20_2_20_points`), by=c("ID"="Sample ID"))

## SU2C cohort----
all( pull( sig.smm.non.null[sig.smm.non.null$cohort%in%c("CTC", "SU2C"),], "ID") %in% clinical.and.terra.ref.pairs.samples$Participant_ID )
# TRUE

sig.smm.non.null <- sig.smm.non.null |> left_join(clinical.and.terra.ref.pairs.samples |> select(Participant_ID, `20 2 20 points`), by=c("ID"="Participant_ID"))

# sig.smm.non.null <- sig.smm.non.null |>
sig.smm.non.null <- sig.smm.non.null |>
  mutate(Mayo_Points=case_when(
    is.na(`20_2_20_points`) ~ `20 2 20 points`,
    is.na(`20 2 20 points`) ~ `20_2_20_points`,
    TRUE ~ -1)) |>
  mutate(IMWG_FISH=Mayo_Points+IMWG_FISH_HighRisk,
         IMWG_Class = case_when(
           IMWG_FISH==0 ~ "Low",
           IMWG_FISH==1 ~ "LowInt",
           IMWG_FISH==2 ~ "Int",
           IMWG_FISH>=3 ~ "High",
           TRUE ~ "NA"
         ))

sig.smm.imwg.stats <- sig.smm.non.null |> group_by(IMWG_Class) |> summarise(n=n(), label=paste0("N=", n), y=9)

sig.smm.non.null |> write_tsv("../data/20241110_MM-like_IMWG_cytogenet_scores.tsv")

table(sig.smm.non.null$IMWG_FISH_HighRisk)

fig.revision.imwg.risk <- ggplot(sig.smm.non.null, aes(factor(IMWG_Class, levels=c("Low", "LowInt", "Int", "High")), total_sig, fill=factor(IMWG_Class, levels=c("Low", "LowInt", "Int", "High")))) + 
  geom_violin(scale = "width", adjust=1.5, alpha=0.25) +
  geom_boxplot(width=0.4) +
  geom_text(data=sig.smm.imwg.stats, aes(y=y, label=label), size=2.5) +
  scale_fill_brewer(palette = 1) +
  scale_y_continuous(breaks=seq(-1, 11)) +
  labs(x="IMWG 2020 Risk\nwith high-risk genomic events", y="MM-like score", fill="") + 
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none")

print(fig.revision.imwg.risk)

ggsave(filename = "../figures/fig.revision.imwg.risk.pdf", plot=fig.revision.imwg.risk, width = 2, height = 2)
ggsave(filename = "../figures/fig.revision.imwg.risk.png", plot=fig.revision.imwg.risk, width = 2, height = 2)
