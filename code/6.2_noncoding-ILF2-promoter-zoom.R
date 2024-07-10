setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

source("0_annotate_samples.R")
# ILF2 hotspots -----------------------------------------------------------

maf.path <- "../data/_MM_Sigs_HDP/20240304_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv.gz"

maf.sbs.2 <- read_tsv(maf.path, num_threads = 7)


ilf2.muts <- maf.sbs.2 |>
  filter(Hugo_Symbol=="ILF2")

positions <- 153643717:153643735
sequence <- strsplit("GTTCTAGATAGCAGAGTGG", split = "")[[1]]

ILF2.REF <- data.frame(Start_position=positions, sequence=sequence)

ilf2.mutants <- ILF2.REF |>
  left_join(ilf2.muts) |>
  filter((!startsWith(Participant_ID, "MMRF") | is.na(Participant_ID))) |> # MMRF will be dealt with later (see 6.3)
  group_by(Start_position, Participant_ID, sub) |>
  slice_head(n=1) |>
  group_by(Start_position, sequence, sub) |>
  summarise(n=sum(!is.na(Tumor_Seq_Allele2)), pts=paste0(Participant_ID, collapse=", ")) |>
  mutate(n=replace_na(n, 0), sub=replace_na(sub, ""))

  
ilf2.mutants.bar <- ilf2.mutants |> 
  ggplot(aes(Start_position, n, fill=sub)) +
  
  geom_bar(stat="identity") +
  geom_text(aes(label=sequence), y=-0.5, size=2.5) +
  
  scale_x_continuous(breaks=c(153643723, 153643730), minor_breaks = 153643717:153643735, labels = scales::comma) +
  scale_y_continuous(breaks=c(0, 3, 6, 9), minor_breaks = 0:10, limits=c(-1,10)) +
  scale_fill_manual(values=c("C>G"="black", "C>T"="red"), na.translate=FALSE) +
  
  labs(x="nucleotide position (chr1, hg19)", y="Mutation counts") +
  
  theme_bw(base_size = 7, base_line_size = 0.25, base_rect_size = 0.25) +
  theme(legend.position = "top",
        panel.grid = element_blank())


ilf2.mutants.bar

ggsave("../figures/ILF2_mutant_plot.png", plot = ilf2.mutants.bar, width = 2, height = 1.5, scale = 1)
ggsave("../figures/ILF2_mutant_plot_legend.pdf", plot = get_legend(ilf2.mutants.bar), width = 2, height = 1.5, scale = 1)
ggsave("../figures/ILF2_mutant_plot_nolegend.pdf", plot = ilf2.mutants.bar +theme(legend.position = "none"), width = 2, height = 1.5, scale = 1)

ilf2.muts.maf <- ilf2.muts |> filter(!startsWith(Participant_ID, "MMRF")) |> filter(Start_position %in% c(153643723, 153643730))

ilf2.muts.maf |> write_tsv("../data/ilf2.muts.maf")


# How many patients -------------------------------------------------------

length(unique(ilf2.muts.maf$Participant_ID))
sort(unique(ilf2.muts.maf$Participant_ID))

# # why this fails?
ilf2.muts.maf |>
  filter(!startsWith(Tumor_Sample_Barcode, "MMRF")) |>
  group_by(Participant_ID, Gender, Age, Stage, HRDTx, Cohort, Assay, IMWG, StageAnd22020) |>
  summarize(Positions=paste(Start_position, collapse = ", ")) |>
  group_by(IMWG) |>
  count()

# how many maf
mm <- read_tsv("../data/export_MM-like_features_MGUS_SMM_MM.tsv") |>
  filter(Cohort%nin%c("MMRF", "Bustoros"))
table(mm$IMWG)

fisher.test((matrix(c(9, 33-9, 0, 177-33), nrow=2)))
# data:  (matrix(c(9, 33 - 9, 0, 177 - 33), nrow = 2))
# p-value = 1.009e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   10.71787      Inf
# sample estimates:
#   odds ratio 
# Inf 


# Clonal ------------------------------------------------------------------

ilf2.muts.maf <- read_tsv("../data/ilf2.muts.maf")

length(unique(ilf2.muts.maf$Tumor_Sample_Barcode))
length(unique(ilf2.muts.maf$Participant_ID))


mutations.ccf.from.su2c <- ilf2.muts.maf |>
  filter(Cohort!="MMRF") |>
  group_by(Start_position, Tumor_Seq_Allele2, Participant_ID, ccf_median) |>
  arrange(Tumor_Sample_Barcode) |>
  slice_head(n=1) |>
  select(Start_position, Tumor_Seq_Allele2, Participant_ID, HRDTx, Stage, ccf_median, SCNA_tau, multiplicity, t_alt_count, t_ref_count) |>
  ungroup()

# stage distribution
mutations.ccf.from.su2c |>
  select(Participant_ID, Stage) |>
  distinct() |>
  group_by(Stage) |>
  count()

# 1q

mutations.ccf.from.su2c

table( mutations.ccf.from.su2c$ccf_median > 0.85 )
table( mutations.ccf.from.su2c$ccf_median > 0.85 )

# CCF ---------------------------------------------------------------------

ilf2.ccf <- ilf2.muts.maf |>
  filter(Start_position %in% c(153643723, 153643730)) |>
  ungroup() |>
  arrange(ccf_CI95_low, ccf_hat) |>
  ggplot(aes(factor(Start_position, levels=c(153643723, 153643730), label=c("153,643,723", "153,643,730")), ccf_hat)) +
  # geom_violin() +
  # geom_jitter(shape=21, width = 0.2, aes(fill=sub)) +
  geom_pointrange(aes(ymin=ccf_CI95_low, ymax=ccf_CI95_high, color=sub), position = position_dodge2(width = 0.9), size=0.25, shape=21) +
  
  # scale_x_discrete(labels=scales::scientific) +
  scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
  scale_color_manual(values=c("C>G"="black", "C>T"="red"), na.translate=FALSE) +
  
  geom_vline(xintercept = 1.5, color="lightgrey") +
  
  labs(x="nucleotide position (chr1, hg19)", y="Cancer cell fraction (%)") +
  
  theme_bw(base_size = 7, base_line_size = 0.25, base_rect_size = 0.25) +
  theme(legend.position = "top",
        panel.grid = element_blank())

ilf2.ccf

ggsave("../figures/ILF2_mutant_plot_ccf.png", plot = ilf2.ccf, width = 2, height = 1.5)
ggsave("../figures/ILF2_mutant_plot_ccf_legend.pdf", plot = get_legend(ilf2.ccf), width = 2, height = 1.5)
ggsave("../figures/ILF2_mutant_plot_ccf_nolegend.pdf", plot = ilf2.ccf+theme(legend.position = "none"), width = 2, height = 1.5)


# Extract patients MUT and NO MUT -----------------------------------------

summary.ilf2 <- ilf2.muts |>
  filter(Start_position %in% c(153643723, 153643730)) |>
  group_by(Participant_ID) |>
  summarise(Start_positions=paste(Start_position, collapse=', '), mutants=paste(sub, collapse=", "))

ilf2.su2c <- clinical.and.terra.ref.pairs.samples |>
  # for ILF 2 mutants
  full_join(summary.ilf2, by=c("HDP_Participant_ID"="Participant_ID")) |>
  filter(Cohort %in% c("CTC", "SU2C")) |>
  # remove MMRF
  filter(!is.na(Participant_ID))  |> 
  select(MRN, PANGEA_ID, CTF_ID, Participant_ID, HDP_Participant_ID, HDP_Reference_Pair, Reference_Pair, `Tissue ID`, Start_positions, mutants )

ilf2.su2c |> write_tsv("../data/ilf2-su2c.ids.tsv")
ilf2.su2c <- read_tsv("../data/ilf2-su2c.ids.tsv")



