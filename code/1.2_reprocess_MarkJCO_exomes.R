setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(readxl)
library(patchwork)

# setwd("code/")

source("0_annotate_samples.R")

### EXPORT ###
# for fig 1
# simple metadata
# mut
# cnv
# sv

# Load Mark Data ----------------------------------------------------------

mark.clinical <- read_xlsx("../../../5_Papers/Internal Manuscripts/Bustoros_JCO_Data/Final Stats for SMM samples_MRNs_Redcap.xlsx")
mark.final.matrix <- read_tsv("../../../5_Papers/Internal Manuscripts/Bustoros_JCO_Data/SMM_study_final_list_Rob_analysis.txt")
mark.seg.file <- read_tsv("../../../5_Papers/Internal Manuscripts/Bustoros_JCO_Data/TN_pairs_95.aggregated.tsv")
mark.maf.file <- read_tsv("../../../5_Papers/Internal Manuscripts/Bustoros_JCO_Data/TN_pairs_95.aggregated_Oncotator_fix_PoN_filter_Bait_Bias_filter_Absolute_Extract.tsv", num_threads = 4)


# Patient filtering -------------------------------------------------------

# 3, 17, 040 due to 
# paste0("SMM_", str_pad( c(3, 17, 40, 74, 9, 29, 37, 38, 42), 3, side = "left", pad= "0"), "_Tumor")
mark.maf.file %>%
  filter(tumor_sample_barcode %in% paste0("SMM-", c("076", "081", "093"), "-Tumor")) %>%
  clipr::write_clip()

# because of doublets with WGS
remove.samples <- paste0("SMM_", c("017", "040"), "_Tumor")
# bring back 003 since WXS is untreated and WGS from same patient is post-EloRd. Removing it from WGS.

# beacuse of no myeloma / IGH mutation or cnv in actual genome data (only relies on FISH from clinic? not good for discovery)
remove.samples <- c(remove.samples, paste0("SMM_", str_pad( c(9, 29, 37, 38, 42), 3,side = "left", pad= "0"), "_Tumor"))
# because of non human contamination
remove.samples <- c(remove.samples, "SMM_074_Tumor")

# additional doublets discovered:
#
# 'SMM-076-Tumor'
# 'SMM-081-Tumor'
# 'SMM-093-Tumor'
# we'll use the baseline WGS from same timepoint

remove.samples <- c(remove.samples, paste0("SMM_", c("076", "081", "093"), "_Tumor"))

mark.final.exomes.matrix <- mark.final.matrix %>%
  filter(str_detect(`Sample ID`, "SMM_[0-9]{3}_Tumor")) %>%
  filter(!(`Sample ID`%in%remove.samples))


# MAF filtering -----------------------------------------------------------

mark.full.maf <- mark.maf.file %>%
  mutate(tumor_sample_barcode=str_replace_all(sample, "-", "_")) %>%
  mutate(Tumor_Sample_Barcode=tumor_sample_barcode, `Sample ID`=Tumor_Sample_Barcode) %>%
  group_by(Tumor_Sample_Barcode, Chromosome, Start_position, End_position, Tumor_Seq_Allele2, sample) %>%
  slice_sample(n=1) # to avoid duplicates

genes.of.interest[genes.of.interest %nin% mark.maf.file$Hugo_Symbol]

# # security check on drivers
# mark.full.maf <- mark.full.maf %>% mutate(Hugo_Symbol = case_when(
#   # Hugo_Symbol=="HIST1H1B"~"H1-5",
#                                  # Hugo_Symbol=="HIST1H1E"~"H1-4",
#                                  # Hugo_Symbol=="RBM16"~"SCAF8",
#                                  # Hugo_Symbol=="KIAA0182"~"GSE1",
#                                  # Hugo_Symbol=="MLL4"~"KMT2D",
#                                  # Hugo_Symbol=="WHSC1"~"NSD2",
#                                  # Hugo_Symbol=="TENT5C"~"FAM46C",
#                                  TRUE~Hugo_Symbol))

mark.final.maf <- mark.full.maf %>%
  filter(tumor_sample_barcode %in% unique(mark.final.exomes.matrix$`Sample ID`))

# duplicated mutations looks clean
mark.final.maf %>% 
  group_by(Chromosome, Start_position) %>%
  summarise(count=n(), pts=paste0(Tumor_Sample_Barcode, collapse = ", ")) %>%
  filter(count>=2) %>%
  arrange(-count)

# this one not clean
filter.non.human.muts <- mark.maf.file %>%
  filter(sample=="SMM-074-Tumor" & Variant_Classification=="Silent") %>%
  select(Chromosome, Start_position)

# map back to other paitnets
COW_MUT_PER_SAMPLE <- filter.non.human.muts %>% 
  left_join(mark.final.maf) %>%
  filter(sample!="SMM-074-Tumor") %>%
  group_by(sample) %>%
  summarise(N_COW=n())

# check non human contamination
mark.maf.file %>%
  group_by(sample) %>%
  summarise(N_SILENT=sum(Variant_Classification=="Silent"),
            N_MISSENSE=sum(Variant_Classification=="Missense_Mutation"),
            RATIO=N_SILENT/N_MISSENSE) %>%
  left_join(COW_MUT_PER_SAMPLE) %>%
  arrange(-N_SILENT)

write_tsv(mark.final.maf, "../data/20230509_bustoros_exome_hg19.maf")

mark.final.maf <- read_tsv("../data/20230509_bustoros_exome_hg19.maf")

mark.final.maf %>%
  select(1:22, Codon_Change, Protein_Change, ref_context, ccf_hat, ccf_CI95_low, ccf_CI95_high) %>%
  write_tsv("../data/20230509_bustoros_exome_light_hg19.maf")

# make seg file look like sample id from paper
mark.final.seg <- mark.seg.file %>%
  mutate(Number = str_match(sample, "SMM-([0-9]{3})_pair")[,2]) %>% 
  mutate(`Sample ID` = paste0("SMM_", Number, "_Tumor")) %>%
  relocate(`Sample ID`) %>%
  select(-c(sample, Number)) %>% 
  filter(`Sample ID`%in%unique(mark.final.exomes.matrix$`Sample ID`))

# export for igv
mark.final.seg %>% write_tsv("../data/20230804_mark_jco_absolute_copy_number.tsv")

mark.final.seg %>%
  mutate(total_cn=rescaled.cn.a1 + rescaled.cn.a2) %>%
  select(1:6, total_cn) %>%
  write_tsv("../data/mark_jco_renorm.seg")


# Make matrix from there --------------------------------------------------

mark.final.exomes.matrix

sort(unique(mark.final.exomes.matrix$`Sample ID`))
sort(unique(mark.final.maf$Tumor_Sample_Barcode))
sort(unique(mark.final.seg$`Sample ID`))

all(sort(unique(mark.final.seg$`Sample ID`))==sort(unique(mark.final.exomes.matrix$`Sample ID`)))
all(sort(unique(mark.final.maf$`Sample ID`))==sort(unique(mark.final.exomes.matrix$`Sample ID`)))


# Structural variants / Matrix report -------------------------------------

# translocations
colnames(mark.final.exomes.matrix)
# t(11;14), t(4;14), t(14;16), t(14;20), t(6;14), Other Tx, MYC translocation

mark.final.exomes.matrix.ready <- mark.final.exomes.matrix |>
  mutate(Gender=SEX,
         Stage="SMM",
         Age=as.character(AGE),
         HRDTx=case_when(`t(4;14)`==1 ~ "t(4;14)",
                         `t(11;14)`==1 ~ "t(11;14)",
                         `t(6;14)`==1 ~ "t(6;14)",
                         `t(14;16)`==1 ~ "t(14;16)",
                         `t(14;20)`==1 ~ "t(14;20)",
                         `HRD (3,5,7,9,11,15,19,21)`==1 ~ "HRD",
                         Other_IgH_translocations==1 ~ "Other tx",
                         TRUE~"Unknown"),
         MYC=case_when( ( `MYC translocation`==1 | `MYC aberrations`==1 ) ~ 1, TRUE ~ 0),
         Cohort="Bustoros",
         Assay="WXS")

# reviewing and seems like
mark.final.exomes.matrix.ready[mark.final.exomes.matrix.ready$`Sample ID`=="SMM_030_Tumor","HRDTx"] <- "HRD"
# others with MRN seem not to have enough cells for FISH or no FISH ordered


## Clinical -----------------------------------------------


clinic.bins <- mark.final.exomes.matrix.ready %>%
  select(`Sample ID`, Gender, Age, Stage, HRDTx, Cohort, Assay) %>%
  pivot_longer(-`Sample ID`)

clinic.bins %>% write_tsv("../data/fig1_tmp/mark_jco_meta.tsv")

clinic.bins %>%
  ggplot(aes(`Sample ID`, name, fill=value)) +
  geom_raster()

genes.of.interest[which(genes.of.interest %nin% mark.final.maf$Hugo_Symbol)]


## Mutation rows -----------------------------------------------------------

genes.of.interest[ genes.of.interest %nin% mark.maf.file$Hugo_Symbol ]

# FIXME investigate mutations that are absent from Mark's dataset (is it due to oncotator?)
all <- expand_grid(`Sample ID`=mark.final.exomes.matrix$`Sample ID`, Hugo_Symbol=genes.of.interest)
mut.bins <- mark.final.maf %>%
  ungroup() %>%
  filter(Hugo_Symbol %in% genes.of.interest & Variant_Classification %in% non.synonymous) %>%
  select(`Sample ID`, Hugo_Symbol, Variant_Classification, Protein_Change, Chromosome, Start_position) %>%
  right_join(all, by=c("Sample ID", "Hugo_Symbol")) %>%
  mutate(name=Hugo_Symbol) %>%
  group_by(`Sample ID`, name) %>% 
  summarise(value = prioritize.variant.classification(Variant_Classification[Variant_Classification!=""]),
            all_mutations = paste0(Variant_Classification[Variant_Classification!=""], collapse = ";", recycle0 = TRUE),
            n_mutations = length(Variant_Classification[Variant_Classification!="" & !is.na(Variant_Classification)])) %>% 
  mutate(value = case_when(is.na(value)|value=="NA"~"",TRUE~value)) %>% 
  mutate(all_mutations = case_when(is.na(all_mutations)|all_mutations=="NA"~"",TRUE~all_mutations)) %>%
  as.data.frame()

# force call from team's annotation / review (2019-2021 exome data)
mut.bins[mut.bins$`Sample ID`=="SMM_019_Tumor" & mut.bins$name=="KRAS", c("value", "all_mutations", "n_mutations")] <- c("Missense_Mutation", "Missense_Mutation", 1)
mut.bins[mut.bins$`Sample ID`=="SMM_093_Tumor" & mut.bins$name=="KRAS", c("value", "all_mutations", "n_mutations")] <- c("Missense_Mutation", "Missense_Mutation", 1)
mut.bins[mut.bins$`Sample ID`=="SMM_056_Tumor" & mut.bins$name=="PRKD2", c("value", "all_mutations", "n_mutations")] <- c("Missense_Mutation", "Missense_Mutation", 1)
mut.bins[mut.bins$`Sample ID`=="SMM_060_Tumor" & mut.bins$name=="PRKD2", c("value", "all_mutations", "n_mutations")] <- c("Missense_Mutation", "Missense_Mutation", 1)
mut.bins[mut.bins$`Sample ID`=="SMM_080_Tumor" & mut.bins$name=="DIS3", c("value", "all_mutations", "n_mutations")] <- c("", "", 0)

mut.bins %>% write_tsv("../data/fig1_tmp/mark_jco_snv.tsv")

mut.bins %>%
  ggplot(aes(`Sample ID`, name, fill=value)) +
  geom_raster()

# mutations (translocations) row

tx.bins <- mark.final.exomes.matrix.ready %>% 
  select(`Sample ID`, `t(11;14)`, `t(14;16)`, `t(14;20)`, `t(4;14)`, `t(6;14)`, `MYC`) %>%
  rename(`t(MYC)`=`MYC`) %>%
  pivot_longer(-`Sample ID`) %>%
  mutate(value=case_when(value==1 ~ "tx",
                         TRUE ~ ""))

tx.bins %>% write_tsv("../data/fig1_tmp/mark_jco_tx.tsv")

cnv.bins <- mark.final.exomes.matrix.ready %>% 
  mutate(Del1p22.1 = `1p_del`,
         Gain1q=`1q_gain`,
         Gain3q26.2=NA, # none in the cohort
         Del6q22.33=`6q_del`,
         Del8p22=del_8p,
         Del13q=`Monosomy 13/ 13q del`,
         Del14q32.32=`14q_del/Monosomy 14`,
         Gain16p13.13=NA, # none in the cohort
         Del17p=`17p_del`,
         Del22q=`del_22q/22`) %>%
  select(`Sample ID`, Del1p22.1, Gain1q, Gain3q26.2, Del6q22.33, Del8p22, Del13q, Del14q32.32, Gain16p13.13, Del17p, Del22q) %>%
  pivot_longer(-`Sample ID`) %>%
  mutate(value=case_when(value==1&str_detect(name, "Gain") ~ "gain",
                         value==1&str_detect(name, "Del") ~ "loss",
                         TRUE ~ ""))

cnv.bins %>% write_tsv("../data/fig1_tmp/mark_jco_cnv.tsv")


# prevalence column 
# will be computed on the fly
prev.bars <- mut.bins %>% 
  group_by(name) %>%
  summarise(Frequency = sum(value %nin% c("", " ", NA, NULL))/n()) %>%
  mutate(pretty.freq = scales::percent(Frequency, accuracy = 0.1))

prev.bars %>%
  ggplot(aes(Frequency, name)) +
  geom_bar(stat="identity", fill="royalblue") +
  geom_text(aes(label=pretty.freq), x=0.5) +
  scale_x_continuous(labels=scales::percent, limits=c(0, NA)) +
  labs(x="Frequency (%)", y="", fill="") +
  theme_bw()
  

# Adding Sept 6 2023 for Sofia clinical -----------------------------------

mark.final.matrix <- read_tsv("../../../5_Papers/Internal Manuscripts/Bustoros_JCO_Data/SMM_study_final_list_Rob_analysis.txt")
mark.final.matrix <- mark.final.matrix |> select(1:5)

# small piece of code for Sofia:
# we want to rescue more genes
more.genes <- c("TRAF3", "CYLD", "BCL7A", "PRDM1", "SETD2", 
           "FGFR3", "IDH1", "HNRNPU", "HNRNPA2B1", "SAMHD1", 
           "IKZF3", "IKBKB", "FIP1L1", "SF3B1", "HIST1H1B", "ACTG1")

genes.to.annotate <- c(genes, all.2)

mark.full.maf <- mark.maf.file %>%
  mutate(tumor_sample_barcode=str_replace_all(sample, "-", "_")) %>%
  mutate(Tumor_Sample_Barcode=tumor_sample_barcode, `Sample ID`=Tumor_Sample_Barcode) %>%
  group_by(Tumor_Sample_Barcode, Chromosome, Start_position, End_position, Tumor_Seq_Allele2, sample) %>%
  slice_sample(n=1)

mark.final.maf <- mark.full.maf %>%
  filter(tumor_sample_barcode %in% unique(mark.final.matrix$`Sample ID`))

at.least.once.more.genes <- more.genes[more.genes %in% mark.maf.file$Hugo_Symbol]
# [1] "PRDM1"  "SETD2"  "SAMHD1" "FIP1L1" "SF3B1"  "ACTG1" 

all <- expand_grid(`Sample ID`=unique(mark.final.matrix$`Sample ID`), Hugo_Symbol=at.least.once.more.genes)

mut.bins <- mark.final.maf %>%
  ungroup() %>%
  filter(Hugo_Symbol %in% at.least.once.more.genes & Variant_Classification %in% non.synonymous) %>%
  select(`Sample ID`, Hugo_Symbol, Variant_Classification, Protein_Change, Chromosome, Start_position) %>%
  right_join(all, by=c("Sample ID", "Hugo_Symbol")) %>%
  mutate(name=Hugo_Symbol) %>%
  group_by(`Sample ID`, name) %>% 
  summarise(value = prioritize.variant.classification(Variant_Classification[Variant_Classification!=""]),
            all_mutations = paste0(Variant_Classification[Variant_Classification!=""], collapse = ";", recycle0 = TRUE),
            n_mutations = length(Variant_Classification[Variant_Classification!="" & !is.na(Variant_Classification)])) %>% 
  mutate(value = case_when(is.na(value)|value=="NA"~"",TRUE~value)) %>% 
  mutate(all_mutations = case_when(is.na(all_mutations)|all_mutations=="NA"~"",TRUE~all_mutations)) %>%
  as.data.frame()