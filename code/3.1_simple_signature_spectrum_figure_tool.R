setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load  -------------------------------------------------------------------

# a simple tool to represent mutation context in patients
# please add to wolF standard pipeline in addition to LEGO plot

library(tidyverse)
library(ComplexHeatmap)
# install.packages("ggpmisc")
library(ggpmisc)

source("0_annotate_samples.R")

HOMEDIR="../data/fig1_tmp/"
METAs=list.files(HOMEDIR, pattern = ".*meta.tsv", full.names = TRUE)
meta.meta <- rbindlist(lapply(METAs, read_tsv))
meta <- meta.meta %>% pivot_wider(id_cols = `Sample ID`)

meta.sigs <- meta %>% filter(Cohort %nin% "Bustoros")

sigs <- read.table("../data/_MM_Sigs_HDP/sig_exposures_final.txt", row.names = 1)
maf.sbs <- tibble(fread("../data/_MM_Sigs_HDP/20231019_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv.gz"), nThread = 4)
# need to remove whatever _1_BM_CD138pos_pair_tumor from the MuTect ID (we already selected 1 sample per patient earlier)


# explore -----------------------------------------------------------------

table(maf.sbs$trinuc_ref, useNA = "ifany") # ref
table(maf.sbs$trinuc_ref_py, useNA = "ifany") # 
table(maf.sbs$trinuc_sub, useNA = "ifany") # substitution with context
table(maf.sbs$sub, useNA = "ifany") # substitution only


sub_ref <- rep( c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each=16)

trinuc_ref_py_ref <- c(
  rep( c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"), 3),
  rep( c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT","GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT"), 3))

skeleton <- tibble(trinuc_ref_py=trinuc_ref_py_ref, sub=sub_ref)

sub.cols <- c("C>A"="dodgerblue", "C>G"="black", "C>T"="red", "T>A"="grey70", "T>C"="olivedrab3", "T>G"="plum2")

maf_to_spectrum <- function(maf) {
  
  maf.to.spectrum <- maf |>
    group_by(trinuc_ref_py, sub) |>
    summarize(n=n()) |>
    right_join(skeleton) |>
    mutate(n=replace_na(n, 0)) |>
    ungroup() |>
    mutate(freq=n/sum(n))
  
  maf.to.spectrum
  
  
}

non.mmrf.maf <- maf.sbs |> filter(!startsWith(Participant_ID, "MMRF"))
non.mmrf.maf.samples <- sort(unique(non.mmrf.maf$Tumor_Sample_Barcode))

pdf("../figures/sbs_plot_su2c.pdf", width=7, height = 2)
for (s in non.mmrf.maf.samples){
  # subset maf
  mini.maf <- non.mmrf.maf |>
    filter(Tumor_Sample_Barcode==s)
  # convert to spectrum
  sp <- maf_to_spectrum(mini.maf)
  #make figure
  pl <- ggplot(sp, aes(trinuc_ref_py, freq, fill=sub)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=sub.cols) +
    facet_grid(~sub, scales = "free_x") +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 1) +
    labs(title = s, x="trinucleotide context", y="frequency") +
    guides(fill = guide_legend(override.aes = list(size = 0.2)))
  # done!
  print(pl)
}

dev.off()

