# code for figure only on the 16 longitudinal samples and MM-like socre
# actual annotation is in Supp Table
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(readxl)

serial.data <- read_xlsx("../DuttaAlberge_CleanClinicalRef.xlsx", sheet="final_serial") # this is now in supplementary tables without RedCap IDs

# Check mutants -----------------------------------------------------------
# (run only once)
source("0_annotate_samples.R")
maf <- read_tsv("../data/_MM_Sigs_HDP/20240304_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv.gz", num_threads = 6)

table(serial.data$patient_name %in% maf$Participant_ID) # all

samples <- sort(unique(maf$Tumor_Sample_Barcode))

tiny.maf <- maf |> 
  filter( Variant_Classification%in%non.synonymous & Hugo_Symbol %in% genes.of.interest & Participant_ID %in% serial.data$patient_name) |>
  select(Participant_ID, Tumor_Sample_Barcode, Tumor_Sample_Barcode_combined, Hugo_Symbol, Chromosome, Start_position, Protein_Change, t_alt_count, t_ref_count, purity, ploidy, ccf_hat)

tiny.maf |> write_tsv("../data/20240313_helper_for_drivers_longitudinal.maf")

# Check SVs ---------------------------------------------------------------

svs <- read_tsv("../data/20230509_all_svs_hg38_liftedOver_fixed_ReferencePair.tsv")
sort(unique(svs$individual))
svs$individual
serial.data$terra_sample_name[!(serial.data$maf_sample_name %in% svs$individual)]
serial.data$terra_sample_name[(serial.data$terra_sample_name %in% svs$individual)]
serial.data$terra_sample_name[(serial.data$maf_sample_name %in% svs$individual)]
serial.data$terra_sample_name[(serial.data$patient_name %in% svs$individual)]

# Reprocess serial data ---------------------------------------------------

# order by time
serial.data <- serial.data |>
  mutate(patient_name = fct_reorder(patient_name, days, max, .desc = TRUE))

# replace whitespace with \n in comments 
serial.data <- serial.data |>
  mutate(comment = str_replace_all(comment, "\\s", "\n"))

# annotate non-stable
serial.data <- serial.data |>
  mutate(isStable = days==0 | comment=="stable")

# rename stable
serial.data <- serial.data |>
  mutate(comment = ifelse(comment=="stable", "s.", comment))

# Figure backbone ---------------------------------------------------------

serial.graph <- 
  ggplot(serial.data, aes(x=days/362.25, y=patient_name)) +
  
  geom_line() +
  geom_point(aes(fill=isStable, shape=sample_type), size=3) +
  
  geom_text(data=subset(serial.data, days==0), aes(label=comment), size=2.5, hjust = 1, position = position_nudge(x=-0.2)) + # initial T0
  geom_text(data=subset(serial.data, days>=1), aes(label=comment), size=2.5, hjust = 0, position = position_nudge(x=0.2)) + # next ones
  
  scale_fill_manual(values=c("TRUE"="black", "FALSE"="#E41A1C")) +
  scale_color_manual(values=c("TRUE"="black", "FALSE"="#E41A1C")) +
  scale_shape_manual(values=c("CTC"=21, "BM"=22)) +
  scale_x_continuous(limits=c(-1, 5), breaks=0:5) +
  
  labs(x="Years since initial biopsy", y="") +
  
  theme_classic() +
  theme(panel.grid = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size=7)) 

serial.graph

ggsave("../figures/serials-graph.pdf", serial.graph, width=3, height = 5)
 