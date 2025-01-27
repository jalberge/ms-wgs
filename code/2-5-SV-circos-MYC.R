# MYC SV summary ----------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

# close neighborhood mutational signature ---------------------------------

sig.neighbors.sv <- read_tsv("../data/SSVGAR/Revisions/data_for_plotting_proportions_of_signatures_subsampled_randomly_for_top_5_genes_for_reviewers (1).tsv")
# sig.neighbors.sv <- read_tsv("../data/SSVGAR/Revisions/data_for_plotting_proportions_of_signatures_FINAL.tsv")
sig.neighbors.sv <- sig.neighbors.sv |> column_to_rownames("...1")

# AID+APOBEC 
100*colSums(sig.neighbors.sv[c( str_c("SBS", c(2, 13, 9, 85)), "N17"),])

  ## reprocess data --------------------------------------------------------

source("0_annotate_samples.R")

final.annot.matrix <- read_tsv("../data/fig1_all_matrix_patients_except_jco_with_disease_stage.tsv")

svs <- read_tsv("../data/20230509_MMRF_SU2C_MARK_OW_SV_filtered_results_baseline_hg38_sd_gt_0.tsv")

length(unique(svs$individual))
# always this painful process of sample to patient matches and across papers
svs <- svs |>
  mutate(Participant_ID=case_when(
    startsWith(individual, "MMRF") ~ str_extract(individual, "MMRF_[0-9]{4}"),
    startsWith(individual, "IID") ~ str_remove(individual, "_T01$"),
    startsWith(individual, "PD") ~ str_remove(individual, "[abc]$"),
    startsWith(individual, "CTF") ~ str_extract(individual, "CTF[0-9]{3}"),
    startsWith(individual, "pM11077BM_138pos19pos") ~ str_extract(individual, "pM11077"),
    startsWith(individual, "pM11272_BMPCs") ~ str_extract(individual, "pM11272"),
                              TRUE ~ individual), .after = individual) |>
  filter(Participant_ID %nin% c("MMRF_1889", "MMRF_2363")) |> # removed from final analysis because inconsistent mutations or copy number
  filter(Participant_ID %nin% c("CSL193")) |> # removed from final matrix because previously in Elo treatment (SMM_TUMOR_003) in JCO Paper
  filter(Participant_ID %nin% c("pM11003")) # removed because re-diagnosed as Amyloidosis

all(svs$Participant_ID %in% final.annot.matrix$Participant_ID)

svs <- svs |> left_join(final.annot.matrix, by="Participant_ID")

myc.low.hg38 <- 127.5E6
myc.high.hg38 <- 128.55E6

#MYC SVs
myc.svs <- svs |> filter(( chr1==8&pos1>myc.low.hg38&pos1<myc.high.hg38 | chr2==8&pos2>myc.low.hg38&pos2<myc.high.hg38) & class=="inter_chr")

ggplot(myc.svs |> filter(chr1==8), aes(pos1)) + geom_vline(xintercept = c(127500000, 129000000)) + geom_histogram()

## classify MYC partners ---------------------------------------
myc.partners <- myc.svs |>
  mutate(gene_MYC=ifelse(chr1==8, gene1, gene2), gene_nonMYC=ifelse(chr1==8, gene2, gene1),
         chr_MYC=ifelse(chr1==8, chr1, chr2), chr_nonMYC=ifelse(chr1==8, chr2, chr1),
         pos_MYC=ifelse(chr1==8, pos1, pos2), pos_nonMYC=ifelse(chr1==8, pos2, pos1)) |>
  mutate(Final_Partner=case_when(
    gene_nonMYC=="GDAP2"~"TENT5C",
    gene_nonMYC=="RPIA"~"IGK",
    gene_nonMYC%in%c("CCR2", "CCR3")~"CCR2/3",
    gene_nonMYC%in%c("FBXW7", "TMEM154", "CENPU")~"FBXW7",
    gene_nonMYC%in%c("BMP6", "TXNDC5", "BLOC1S5", "MUTED", "EEF1E1", "SLC35B3")~"TXNDC5",
    gene_nonMYC%in%c("FOXO3", "ARMC2", "SESN1")~"FOXO3",
    gene_nonMYC%in%c("FCHSD2", "P2RY2")~"P2RY2",
    gene_nonMYC%in%c("TEX22","C14orf80","TEDC1","TMEM121","KIAA0125")~"IGH",
    gene_nonMYC%in%c("PRAME","GGTLC2","MIR650","IGLL5","RSPH14")~"IGL",
                   TRUE~gene_nonMYC)) |>
  group_by(Final_Partner, Participant_ID) |>
  slice_head(n=1) |>
  group_by(Final_Partner) |>
  mutate(Count_Participants=length(unique(Participant_ID))) |>
# Problem with particiapnts with multiple Final Partner....
  group_by(Participant_ID, Stage) |>
  summarize(Final_Partner=paste(Final_Partner, collapse = ", "))

## compute IG vs MM partners ------------------------------------------

myc.partners.ig.annot <- myc.partners |>
  # filter(Count_Participants >= 3) |>
  mutate(PartnerClass=case_when(
    str_detect(Final_Partner, "IG[HKL]") ~ "Ig",
    TRUE ~ "Non-Ig")) |>
  mutate(Side = case_when(Stage %in% c("MGUS", "SMM") ~ -0.01, TRUE ~ -1))

## are there more MYC in SMM than in MGUS ------------------

# total MGUS genomes=37 (with oben ad wtc)
# total SMM genomes=120 (same)
ig.vs.non.ig.myc.tx <- myc.partners.ig.annot |>
  mutate(isMM=Stage=="MM") |>
  group_by(isMM, PartnerClass) |>
  count() |>
  pivot_wider(id_cols = PartnerClass, names_from = isMM, values_from = n) |>
  tibble::column_to_rownames("PartnerClass")
# NUMBERS IN THE PAPER
ig.vs.non.ig.myc.tx

fisher.test( ig.vs.non.ig.myc.tx )
# exactly comparable

# total.myc.prop
# this is total number of MYC in MGUS/SMM
# 138 MGUS/SMM from PCROWD + 4 WTC + 15 Oben
prop.test(9, 157)
# 1-sample proportions test with continuity correction
# 
# data:  9 out of 157, null probability 0.5
# X-squared = 121.3, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5
# 95 percent confidence interval:
#   0.02823141 0.10931946
# sample estimates:
#   p 
# 0.05732484 
## are there more IG MYC in MM than MGUS/SMM? ------------------

# total MGUS genomes=37 (with oben ad wtc)
# total SMM genomes=120 (same)
nrow(myc.partners.ig.annot)
length(unique(myc.partners.ig.annot$Participant_ID))

# total MYC MGUS = 1/37
# total MYC SMM = 8/120
# total MYC MM = 151/812
fisher.test(matrix(c(1, 37-1, 8, 120-8, 151, 812-151), nrow = 2))
fisher.test(matrix(c(1, 37-1, 8, 120-8), nrow = 2))
fisher.test(matrix(c(8, 120-8, 151, 812-151), nrow = 2))

mgus.vs.smm.myc.tx <- myc.partners.ig.annot |>
  filter(Stage != "MM") |>
  group_by(Stage, PartnerClass) |>
  count() |>
  pivot_wider(id_cols = PartnerClass, names_from = Stage, values_from = n) |>
  tibble::column_to_rownames("PartnerClass") |>
  mutate(across(everything(), ~ replace_na(., 0)))
# NUMBERS IN THE PAPER
mgus.vs.smm.myc.tx

fisher.test( mgus.vs.smm.myc.tx )
# exactly comparable

# total.myc.prop
# this is total number of MYC in MGUS/SMM
# 138 MGUS/SMM from PCROWD + 4 WTC + 15 Oben
prop.test(9, 157)
# 1-sample proportions test with continuity correction
# 
# data:  9 out of 157, null probability 0.5
# X-squared = 121.3, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5
# 95 percent confidence interval:
#   0.02823141 0.10931946
# sample estimates:
#   p 
# 0.05732484 


# how many gain8q in MGUS/SMM pts -----------------------------------------

mm <- read_tsv("../data/export_MM-like_features_MGUS_SMM_MM.tsv")
#   df$amp_chr_8q <- ifelse(df$del_G.peak.17_del_8q24.21 >=1 | 
#                           df$del_G.peak.18_del_8q24.21 >=1 | 
#                           df$amp_G.peak.3_amp_8q24.21 >=1 | 
#                           df$MYC_translocation >=1, 1, 0)

mm$MYC_FOCAL = 1 * ( mm$del_G.peak.17_del_8q24.21 >=1 | mm$del_G.peak.18_del_8q24.21 >=1 | mm$amp_G.peak.3_amp_8q24.21 >=1 )
mm$MYC_translocation = 1 * (mm$MYC_translocation >= 1)

mm.with.myc_tx <- mm |>
  filter(Cohort %in% c("CTC", "MMRF", "OBEN", "SU2C", "WTC")) |>
  left_join(myc.partners.ig.annot |> mutate(MYC_TX=1), by=c("Sample ID"="Participant_ID"))

mm.with.myc_tx$largeTX = 1 * (mm.with.myc_tx$MYC_translocation >= 1 | mm.with.myc_tx$MYC_TX >= 1)

# 7/9 MGUS/SMM with translocation have a focal gain/loss

# ID	Sex	Age	Disease_Status	cohort	total_sig	amp_chr_8q	MYC_TX	Comment
# CSL201	F	34	SMM	SU2C	1	0	1	In Breakpointer; extremely focal gain; is it missing from the survival data?
#   CTF023	M	84	SMM	CTC	3	0	1	In Breakpointer; extremely focal gain; is it missing from the survival data?
#   CTF038	M	62	SMM	CTC	3	0	1	In Breakpointer; extremely focal gain; is it missing from the survival data?
#   CTF025	F	77	SMM	CTC	4	1	1	Is it focal? Yes
# DL2029	M	54	SMM	SU2C	4	1	1	Is it focal? No
# IL0291	M	70	SMM	SU2C	4	1	1	Is it focal? Yes
# PD47570	M	50	MGUS	OBEN	2	1	1	Is it focal? Yes
# pM10109	M	47	SMM	SU2C	5	1	1	Is it focal? No
# pM9660	F	55	SMM	SU2C	4	1	1	Is it focal? Yes
mm.like <- read_tsv("../data/export_MM-like_score.tsv")

mgus.smm.genomes.score <- mm.like |> 
  filter(Disease_Status %in% c("MGUS", "SMM") & 
           cohort %in% c("CTC", "MMRF", "OBEN", "SU2C", "WTC"))
# the 157 MGUS/SMM WGS stays consistent here. (checked)


# FIXME -------------------------------------------------------------------

# left_join fails (to investigate)
discrepant.tx.cna <- mgus.smm.genomes.score |>
  left_join(myc.partners.ig.annot |> mutate(MYC_TX=1), by=c("ID"="Participant_ID")) |>
  mutate(MYC_TX=replace_na(MYC_TX, 0))
  # filter(MYC_TX != amp_chr_8q)
# 41 MM with MYC Tx but Without Gain/Loss not counted in the score? Double checked
# are these local SV, without transloctaion?
discrepant.tx.cna |> select(1:5, total_sig, amp_chr_8q, MYC_TX, Final_Partner, PartnerClass ) |> write_tsv("../data/discrepant_MYC_Tx_to_review.tsv")




stage.wgs.counts <- final.annot.matrix |> filter(!startsWith(Participant_ID, "SMM")) |> group_by(Stage) |> count()

myc.stage.numbers <- myc.partners.ig.annot |>
  group_by(Stage) |>
  summarise(nn=n()) |>
  left_join(stage.wgs.counts, by=c("Stage"))

myc.stage.numbers
# > myc.stage.numbers
# # A tibble: 3 × 3
# Stage    nn     n
# <chr> <int> <int>
#   1 MGUS      1    37
# 2 MM      151   812
# 3 SMM       8   120
myc.stage.binom.confit <- myc.stage.numbers |>
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) |>
  group_by(Stage) |>
  do(tidy(prop.test(x = .$nn, n = .$n, p=0.8, alternative = "two.sided")))
myc.stage.binom.confit

myc.stage.numbers |> 
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) |>
  arrange(Stage) |>
  column_to_rownames("Stage") |>
  DescTools::CochranArmitageTest()

initial.disease.stage.palette <- c("MM"="#253494", 
                                   "SMM\nProgressed"="#225ea8", 
                                   "SMM"="#41b6c4", 
                                   "MGUS"="#7fcdbb")

myc.stage.binom.confit |> write_tsv("../data/myc_counts_confint.tsv")

## percent MYC barplot -----------------------------------------------------

myc.stage.binom.confit <-  read_tsv("../data/myc_counts_confint.tsv")

initial.disease.stage.palette <- c("MM"="#253494", 
                                   "SMM\nProgressed"="#225ea8", 
                                   "SMM"="#41b6c4", 
                                   "MGUS"="#7fcdbb")

percent.myc <- ggplot(myc.stage.binom.confit, aes(Stage, estimate, fill=Stage)) +
  
  geom_bar(stat="identity") +
  geom_linerange(aes(ymin=conf.low, ymax=conf.high)) +
  
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=initial.disease.stage.palette) +
  
  labs(x="", y="Fraction positive for MYC") +
  
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        legend.position="none")

percent.myc

ggsave("../figures/percent.myc.pdf", width = 2, height = 2)
ggsave("../figures/percent.myc.png", width = 2, height = 2)

## show lollipop plot ------------------------------------------

# progressors ?
myc.partners.ig.annot |> 
  filter(Stage%in% c("MGUS", "SMM")) |> 
  pull(Participant_ID) |> 
  sort() |> 
  dput()

progressions <- c("CSL201"="NP", 
  "CTF023"="NP", 
  "CTF025"="NP", 
  "CTF038"="NP", 
  "DL2029"="NP", 
  "IL0291"="P", 
  "PD47570"="P", 
  "pM10109"="NP", 
  "pM9660"="NP") |>
  as.list() |>
  data.frame() |>
  pivot_longer(cols=everything(), names_to = "Participant_ID", values_to = "Progression")
  
myc.partners.ig.annot <- myc.partners.ig.annot |> left_join(progressions)

myc.partners.ig.annot |> 
  filter(!is.na(Progression)) |>
  group_by(Participant_ID, Progression, Final_Partner) |>
  count() |>
  View()

myc.lollipop <- ggplot(myc.partners.ig.annot |> filter(Stage %in% c("MGUS","SMM")), aes(x = pos_MYC)) + 
  stat_bin(aes(y=..count../sum(..count..)),geom="step", data=myc.partners.ig.annot |> filter(Stage=="MM"), color="#888888") +
  geom_segment(aes(xend=pos_MYC, y=0, yend=Side)) + 
  geom_point(aes(y=Side, color=Progression, fill=Progression), shape=21) +
  geom_text(aes(label=Final_Partner, y=-0.015), angle=90, hjust=1, fontface="italic") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits=c(-0.05, 0.12)) +
  scale_fill_manual(values=c("P"="black", "NP"="white")) +
  scale_color_manual(values=c("P"="black", "NP"="black")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top")

myc.lollipop
ggsave("../figures/myc.lollipop.png", myc.lollipop, width = 2, height = 1, scale = 3.5)
ggsave("../figures/myc.lollipop.pdf", myc.lollipop, width = 2, height = 1, scale = 3.5)
  
## summarize at partner level  ---------------------------------

myc.summary <- myc.partners |>
  group_by(Final_Partner) |>
  summarise(Participant_IDs=paste(Participant_ID, collapse=","), 
            Chr=first(chr_nonMYC),
            Coordinate=median(pos_nonMYC),
            Count_SV=n(),
            Count_Participants=length(unique(Participant_ID)),
            Count_MGUS=sum(Stage=="MGUS"),
            Count_SMM=sum(Stage=="SMM"),
            Count_MM=sum(Stage=="MM")
            )



## transform for circos --------------------------------------------------


# notes from Barwick Light Chain
# MYC translocations to non-immunoglobulin genes 
# [t(other); N = 81], 
#. IgHMYC translocations [t(8;14); N= 27], 
#. IgK-MYC translocations [t(2;8); N= 13], 
#. and IgL-MYC translocations [t(8;22); N = 25].

set.seed(123)
myc.summary <- myc.summary |> 
  arrange(Count_Participants)
  # filter(Count_Participants >= 3)
bed1 <- myc.summary |> transmute(chr=paste0("chr", 8), start=127733434)#, end=127744951)
bed2 <- myc.summary |> transmute(chr=paste0("chr", Chr), start=Coordinate)#, end=Coordinate+500000)
bed.text <- myc.summary |> transmute(chr=paste0("chr", Chr), 
                                     start=Coordinate, 
                                     end=Coordinate, 
                                     value=1, 
                                     name=case_when(Count_Participants>=3~ Final_Partner, TRUE ~ ""), 
                                     color=case_when(Final_Partner%in% c("IGH", "IGK", "IGL")~"#111111", TRUE ~ "#AAAAAA"))
link.line.width <- myc.summary |> 
  mutate(LineWidth=Count_Participants/10) |>
  pull(LineWidth)#, end=Coordinate+500000)
link.color <- myc.summary |> 
  mutate(color=case_when(Count_Participants >= 3 ~ "salmon", TRUE ~ "#AAAAAA")) |>
  pull(color)


## draw circos plot  -----------------------------------------------------------


# adj=c(1, 0.5)
pdf("../figures/myc-mgus-smm-mm-circosplot.pdf", width = 4, height = 4)
circos.clear()
set.seed(123)
circos.par("start.degree" = 90, "circle.margin"=0.3)
circos.initializeWithIdeogram(plotType = c("ideogram", "labels"), 
                              species="hg38", 
                              chromosome.index = paste0("chr", 1:22), 
                              labels.cex = 0.7)
circos.genomicLink(bed1, bed2, lwd=link.line.width, col = link.color)
circos.genomicTrack(bed.text, track.index = get.current.track.index(), ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 1 + mm_y(10), cex = 0.7, labels.column = 2, col = value$color, facing = "clockwise", niceFacing = TRUE, ...)
}, bg.border = NA)
dev.off()

