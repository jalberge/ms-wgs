# this contains :
# - gathering q values for SNV/indels; CNV; SV drivers
# - testing for enrichment of drivers in MM vs MGUS/SMM
# - creating the CoMut plot of the 1,030 individuals
# noting raw data is not publicly available (contains RedCapIDs and list of somatic mutations) - will make it controlled-access

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(data.table)
library(RColorBrewer)

source("0_annotate_samples.R")
HOMEDIR="../data/fig1_tmp/"
# input data columns must be `Sample ID`, name, value

# reclassify hyperdiploids based on andrea --------------------------------

## Part 1 hrd number of chromosomes
## reannotate hyperdiploidy based on overall Andrea's re-analysis of CNV based on total copy number ratios
cnv.new <- read_tsv("../Andrea/SU2C_CoMM/results/Paper_analysis/SU2C-CoMM-Bus_1136samp_complete_CN_callset_102323.txt")

# rename samples from Andrea
cnv.new <- cnv.new %>%    
  mutate(sample = case_when(
    str_starts(sample, "CTF") ~ str_extract(sample, "CTF[0-9]{1,}"),
    str_starts(sample, "IID") ~ str_extract(sample, "IID_H19606[1-4]"),
    str_starts(sample, "MMRF") ~ str_extract(sample, "MMRF_[0-9]{4}"),
    str_starts(sample, "PD") ~ str_extract(sample, "PD[0-9]{5}"),
    str_starts(sample, "MBp") ~ str_extract(sample, "MBp[0-9]{2}"),
    str_starts(sample, "SMM-") ~ paste("SMM", str_extract(sample, "[0-9]{1,}"), "Tumor", sep = "_"),
    str_ends(sample, "BM_138pos19pos") ~ str_remove(sample, "BM_138pos19pos"), # quick hack to remove RedCap IDs before going public
    str_starts(sample, "pM") & str_ends(sample, "_BMPCs") ~ str_extract(sample, "^pM[0-9]{5}"),  # quick hack to remove RedCap IDs before going public
    TRUE ~ sample))

jb.hrd.calls <- hrd.jb(cnv.new)
hrd.counts <- jb.hrd.calls %>% 
  select(sample, total_HD_chr) %>% 
  pivot_wider(names_from = sample, values_from = total_HD_chr) %>%
  as.matrix()
rownames(hrd.counts) <- "HRD_chr_count"
hrds <- names(hrd.counts[,hrd.counts>=2])

hrd.counts %>% as.data.frame %>% tibble::rownames_to_column("Variable")%>% write_tsv("../data/FINAL_Nov2_HRD_Counts_Numeric_Matrix.tsv")


# META --------------------------------------------------------------------

# CLEAN and SAVE this
meta <- read_tsv("../data/FINAL_meta_keys.tsv")
meta.df <- meta %>% 
  mutate(
    IMWG=case_when(
      HRDTx %in% c("t(11;14)", "t(6;14)", "t(12;14)") ~ "Cyclin D",
      HRDTx %in% c("t(4;14)") ~ "MMSET",
      HRDTx %in% c("t(14;16)", "t(14;20)", "t(8;14)") ~ "MAF",
      HRDTx == "HRD" & `Sample ID` %nin% hrds ~ "Unclassified", # 1 case from JCO Mark paper where our CN analysis doesn't find hyperdiploidy
      HRDTx %in% c("HRD") | `Sample ID` %in% hrds ~ "Hyperdiploid", # few cases where andrea can rescue hyperdiploid at low purity
      TRUE ~ "Unclassified"),
    StageAnd22020=Stage,
    Stage=case_when(Stage %in% c("LRSMM", "IRSMM", "HRSMM")~"SMM", str_detect(Stage, "NDMM")~"MM", TRUE ~ Stage),
    Cohort=case_when(Cohort=="MARK"~"SU2C", TRUE~Cohort)) %>%
  mutate(IMWG=factor(IMWG, levels=c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) %>%
  mutate(Gender=factor(Gender, levels=c("F", "M"))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  tibble::column_to_rownames("Sample ID") %>%
  as.data.frame()

meta.df |> tibble::rownames_to_column("Sample ID") %>% write_tsv("../data/FINAL_meta_columns_matrix.tsv")


# P values ----------------------------------------------------------------

## for initiating tx ------------------------------------------------------
tx.order <- c("Trisomies", "t(11;14)(CCND1)", "t(6;14)(CCND3)", "t(12;14)(CCND2)", "t(4;14)(MMSET)", "t(14;16)(MAF)", "t(14;20)(MAFB)", "t(8;14)(MAFA)", "t(MYC;IgH/K/L)")
q.tx <- data.frame(row.names = tx.order,
                   q=rep(x=1, length(tx.order)),
                   group=rep(x="Translocation", length(tx.order)))

## for snvs ---------------------------------------------------------------

# export supp table 2
# bind_rows(genes.of.interest.table, genes.of.interest.rare.table) |>
  # arrange(q) |>
  # clipr::write_clip()

# how many?
nrow(genes.of.interest.rare.table)+nrow(genes.of.interest.table) -1 # IGLL5
# 61

# how many rare 
nrow(genes.of.interest.rare.table)

# then filter with higher prevalence used in matrix
# genes of interest only >1%
q.snv <- genes.of.interest.table %>% 
  filter(gene!="IGLL5") %>%
  select(gene, q) %>% 
  tibble::column_to_rownames("gene") %>% 
  mutate(group="Point mutation\nSmall indel") %>% 
  as.data.frame() 

# rare are in genes.of.interest.rare.table
# View(genes.of.interest.rare.table)

## for cnvs ---------------------------------------------------------------

# first extract significant values that are prevalent enough
cnv.threshold <- 0.03

list.of.cnv.above.threshold.frequency <- cnv.new %>%
  summarise(across(where(is.numeric), ~ sum(. > 0.14) / nrow(cnv.new))) %>% # corresponds to subclonal detection 0.1 in log2 ratio space
  pivot_longer(cols = everything()) %>%
  filter(value>=cnv.threshold) %>%
  pull(name)

# move to git-hosted GISTIC files
# gistic.broad <- read_tsv("../Andrea/SU2C_CoMM/results/GISTIC_runs_results/539216_merge_CoMM_SU2C_q10/CoMM-SU2C_q10.broad_significance_results.txt")
gistic.broad <- read_tsv("../annot/gistic2_broad_results.txt")
# gistic.focal <- read_tsv("../Andrea/SU2C_CoMM/results/GISTIC_runs_results/539216_merge_CoMM_SU2C_q10/CoMM-SU2C_q10.all_lesions.conf_95.txt")
gistic.focal <- read_tsv("../annot/gistic2.all_lesions.conf_95.txt")

# don't need to report hrd arms / avoid double-counting the same intiaiting event
broad.significant <- c("amp_chr_1q", "amp_chr_6p", "amp_chr_18p", "amp_chr_18q", "del_chr_8p", "del_chr_13q", "del_chr_16q")

broad.cn.q.values <- gistic.broad %>%
  select(Arm, `Amp q-value`, `Del q-value`) %>%
  pivot_longer(cols = -Arm) %>%
  mutate(cn_name = case_when(
    startsWith(name, "Amp")~paste0("amp_chr_", Arm),
    startsWith(name, "Del")~paste0("del_chr_", Arm),
    TRUE ~ "FAIL")) %>%
  filter(cn_name %in% broad.significant) %>%
  transmute(q=value, cn_name)

focal.cn.q.values <- gistic.focal %>%
  filter(!endsWith(`Unique Name`, "CN values")) %>%
  select(`Unique Name`, Descriptor, `q values`) %>%
  mutate(peak_number=str_extract(`Unique Name`, "[0-9]{1,}")) %>%
  mutate(andrea_peak_number=row_number()) %>%
  mutate(cn_name=case_when(
    startsWith(`Unique Name`, "Amplification") ~ paste0("amp_G-peak-", andrea_peak_number, "_amp_", Descriptor),
    startsWith(`Unique Name`, "Deletion") ~ paste0("del_G-peak-", andrea_peak_number, "_del_", Descriptor),
    TRUE ~ "FAIL"
  )) %>%
  filter(cn_name %in% list.of.cnv.above.threshold.frequency) %>%
  transmute(q=`q values`, cn_name)

q.cnvs.names <- rbind(focal.cn.q.values, broad.cn.q.values) %>%
  arrange(q)

# rename cnv
q.cnvs.names <- q.cnvs.names %>%
  mutate(cn_name=str_replace(cn_name, "del_G-peak-17_del_8q24.21", "del_G-peak-17_del_8q24.21a")) %>%
  mutate(cn_name=str_replace(cn_name, "del_G-peak-18_del_8q24.21", "del_G-peak-18_del_8q24.21b")) %>%
  mutate(short_name=str_replace(cn_name, "amp_chr_", "+")) %>%
  mutate(short_name=str_replace(cn_name, "amp_chr_", "+")) %>%
  mutate(short_name=str_replace(short_name, "del_chr_", "-")) %>%
  mutate(short_name=str_replace(short_name, "del_G-peak-[0-9]{1,}_del_", "-")) %>%
  mutate(short_name=str_replace(short_name, "amp_G-peak-[1234]_amp_", "+"))


q.cnvs <- q.cnvs.names %>%
  select(q, short_name) %>%
  tibble::column_to_rownames("short_name") %>% 
  mutate(group="Copy number\nAbnormality") %>% 
  as.data.frame()


## for svs ----------------------------------------------------------------

fisher.method <- function(p) pchisq(-2 * sum(log(p)), 2*length(p), lower.tail = FALSE)

# move to final annot
# q.svs <- read_tsv("../data/Xavi SV LOF results.tsv") %>% transmute(gene, q=q_value_LOF) %>% filter(q<0.1)
q.svs <- read_tsv("../annot/SSVGAR_v11_hits.tsv") %>% transmute(gene, q=q_value_LOF) %>% filter(q<0.1)

q.svs <- q.svs %>%
  mutate(gene=ifelse(gene=="RB1", "RB1 SV", gene)) %>% # to avoid mismatch with RB1 from MutSig2CV or GISTIC2
  mutate(gene=ifelse(gene=="WWOX", "WWOX (MAF)", gene)) %>% # this is actually MAF translocated patients - WWOX is a known fragile site. Keep both annotations.
  add_row(gene = "SP140 / SP140L", q = fisher.method(q.svs %>% filter(gene %in% c("SP140", "SP140L")) %>% pull(q))) %>% # clustered events
  add_row(gene = "TRAF3 / AMN", q = fisher.method(q.svs %>% filter(gene %in% c("TRAF3", "AMN")) %>% pull(q))) %>%  # clustered events
  add_row(gene = "PRSS2 / MTRNR2L6", q = fisher.method(q.svs %>% filter(gene %in% c("PRSS2", "MTRNR2L6")) %>% pull(q))) %>%  # clustered events
  arrange(q) %>%
  tibble::column_to_rownames("gene") %>%
  mutate(group="Structural\nVariant") %>%
  as.data.frame()

## group q val ------------------------------------------------------------

gene.q.val <- rbind(q.tx, q.snv, q.cnvs, q.svs)
gene.q.val$group <- factor(gene.q.val$group, levels=c("Translocation", "Point mutation\nSmall indel", "Copy number\nAbnormality", "Structural\nVariant"))

# SNV ---------------------------------------------------------------------

# give an arbitrary score to mutations based on the variant_type annotation
# this is only useful if clustering individuals (same as DLBCL paper). 
# We may not use it here. (standard CoMut plot instead)

SNVs=list.files(HOMEDIR, pattern = ".*snv.tsv", full.names = TRUE)
final.snv <- rbindlist(lapply(SNVs, read_tsv))
final.snv <- final.snv %>% filter(name!="IGLL5") # exon 2-3 are variable IGL region
final.snv.numeric.matrix <- final.snv  %>% 
  mutate(value = case_when(
    n_mutations > 1 ~ 3, # should it be the sum of unique scores?
    value %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", # protein truncating get an arbitrary score of 3 
                 "Frame_Shift_Del", "Start_Codon_Del", "Start_Codon_SNP", 
                 "START_CODON_SNP", "Stop_Codon_Del") ~ 3,
    value %in% c("COULD_NOT_DETERMINE", "DE_NOVO_START_IN_FRAME", "Frame_Shift_Ins",
                 "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", 
                 "Nonstop_Mutation", "Splice_Site", "Stop_Codon_Del") ~ 2,
    TRUE ~ 0 # Silent mutations have been removed ~ should we bring them back?
  )) %>%
  pivot_wider(id_cols = `Sample ID`, values_from = value, names_from = name) %>%
  tibble::column_to_rownames("Sample ID") %>%
  as.matrix() %>%
  t()

final.snv.numeric.matrix <- final.snv.numeric.matrix[rownames(q.snv),] # reorder rows
final.snv.numeric.matrix %>% as.data.frame %>% tibble::rownames_to_column("Variable") %>% write_tsv("../data/FINAL_Oct9_SNV_Hotspot_Numeric_Matrix.tsv")

# SV ----------------------------------------------------------------------

# give an arbitrary score to mutations based on the variant_type annotation
# this is only useful if clustering individuals (same as DLBCL paper). 
# We may not use it here. (standard CoMut plot instead)

final.sv <- read_tsv("../data/FINAL_Jun17_SV_Hotspot_Matrix.tsv")
# final.sv <- read_tsv("../data/FINAL_Oct9_SV_Hotspot_Matrix.tsv")

final.sv.numeric.matrix <- final.sv %>%
  mutate(value=case_when(
    value=="At least one SV present; at least one LoF" ~ 2,
    value=="At least one SV present; no LoF" ~ 1,
    TRUE ~ 0
  )) %>%
  pivot_wider(id_cols = `Sample ID`, values_from = value, names_from = name) %>%
  tibble::column_to_rownames("Sample ID") %>%
  as.matrix() %>%
  t()

rownames(final.sv.numeric.matrix)[rownames(final.sv.numeric.matrix) == "WWOX"] <- "WWOX (MAF)"
final.sv.numeric.matrix <- final.sv.numeric.matrix[rownames(q.svs)[rownames(q.svs) %in% rownames(final.sv.numeric.matrix)],] # reorder rows
final.sv.numeric.matrix %>% as.data.frame %>%  tibble::rownames_to_column("Variable") %>% write_tsv("../data/FINAL_Oct9_SV_Hotspot_Numeric_Matrix.tsv")

# CNV ---------------------------------------------------------------------

# give an arbitrary score to mutations based on the variant_type annotation
# this is only useful if clustering individuals (same as DLBCL paper). 
# We may not use it here. (standard CoMut plot instead)

# to check
# final.cnv.numeric.matrix <- cnv %>%
final.cnv.numeric.matrix <- cnv.new %>%
  select(1, any_of(intersect( c(focal.cn.q.values$cn_name, broad.cn.q.values$cn_name), list.of.cnv.above.threshold.frequency))) %>%
  tibble::column_to_rownames("sample") %>%
  mutate_all(~ case_when(. >= 2 ~ 3,. >= 1 ~ 2, . >= 0.14 & . < 1 ~ 1, TRUE ~ 0)) %>%
  as.matrix() %>%
  t()

# two focal dels around MYC on the same cytoband - renaming a and b
rownames(final.cnv.numeric.matrix)[which( rownames(final.cnv.numeric.matrix) == "del_G-peak-17_del_8q24.21" ) ] <- "del_G-peak-17_del_8q24.21a"
rownames(final.cnv.numeric.matrix)[which( rownames(final.cnv.numeric.matrix) == "del_G-peak-18_del_8q24.21" ) ] <- "del_G-peak-18_del_8q24.21b"

switch.names <- q.cnvs.names %>% tibble::column_to_rownames("cn_name")
rownames(final.cnv.numeric.matrix) <- switch.names[ rownames(final.cnv.numeric.matrix), "short_name"] 

final.cnv.numeric.matrix <- final.cnv.numeric.matrix[rownames(q.cnvs %>% arrange(q)),] # reorder rows
final.cnv.numeric.matrix %>% as.data.frame %>% tibble::rownames_to_column("Variable") %>% write_tsv("../data/FINAL_Nov2_CNV_Focal_Numeric_Matrix.tsv")

# Tx subgroup -------------------------------------------------------------

# give an arbitrary score to mutations based on the variant_type annotation
# this is only useful if clustering individuals (same as DLBCL paper). 
# We may not use it here. (standard CoMut plot instead)

TXs=list.files(HOMEDIR, pattern = ".*tx.tsv", full.names = TRUE)
tx.tx <- rbindlist(lapply(TXs, read_tsv))
final.tx.numeric <- tx.tx %>% 
  mutate(value=case_when(
    value=="tx"&!is.na(value)&name=="t(MYC)"~5, # discuss how to adjust this number
    value=="tx"&!is.na(value)&name!="t(MYC)"~5,
    is.na(value)~0,
    TRUE ~ 0
  )) %>%
  mutate(name=case_when(
    name=="t(11;14)"~"t(11;14)(CCND1)",
    name=="t(14;16)"~"t(14;16)(MAF)",
    name=="t(14;20)"~"t(14;20)(MAFB)",
    name=="t(4;14)"~"t(4;14)(MMSET)",
    name=="t(6;14)"~"t(6;14)(CCND3)",
    name=="t(MYC)"~"t(MYC;IgH/K/L)",
    name=="t(12;14)"~"t(12;14)(CCND2)",
    name=="t(8;14)"~"t(8;14)(MAFA)",
    TRUE~name)) %>%
  pivot_wider(id_cols = `Sample ID`) %>%
  tibble::column_to_rownames("Sample ID") %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  as.matrix() %>%
  t()

final.tx.numeric %>% as.data.frame %>% tibble::rownames_to_column("Variable")%>% write_tsv("../data/FINAL_Jan4_TX_Focal_Numeric_Matrix.tsv")


# Upset to check patient overlap ------------------------------------------

lt <- list(
  HRD = colnames(hrd.counts),
  CNV = colnames(final.cnv.numeric.matrix),
  SNV = colnames(final.snv.numeric.matrix),
  SV = colnames(final.sv.numeric.matrix),
  TX = colnames(final.tx.numeric))
comparison.lt <- list_to_matrix(lt)

comparison.lt %>% as.data.frame() %>% filter(SNV==1 & SV==1 & TX==1 & CNV ==0) %>% rownames()
comparison.lt %>% as.data.frame() %>% filter(SNV==0 & SV==0 & TX==0 & CNV ==1) %>% rownames()

matrix.samples <- make_comb_mat(comparison.lt)
UpSet(matrix.samples)

# final.clean.samples <- comparison.lt %>% as.data.frame() %>% filter(SNV==1 & SV==1 & CNV == 1 & TX ==1) %>% rownames()
final.clean.samples <- comparison.lt %>% as.data.frame() %>% filter(SNV==1 & CNV == 1 & TX ==1) %>% rownames()

# Plot --------------------------------------------------------------------

# empty SV - ideally set to mean of each feature?
nrow.sv <- nrow(final.sv.numeric.matrix)
ncol.sv <- sum(!final.clean.samples %in% colnames(final.sv.numeric.matrix))
add.empty.sv <- matrix(0, nrow=nrow.sv, ncol = ncol.sv)
rownames(add.empty.sv) <- rownames(final.sv.numeric.matrix)
colnames(add.empty.sv) <- final.clean.samples[!(final.clean.samples %in% colnames(final.sv.numeric.matrix))]
final.sv.numeric.matrix <- cbind(final.sv.numeric.matrix, add.empty.sv)[rownames(final.sv.numeric.matrix), final.clean.samples]

# Grouped together --------------------------------------------------------

final.numeric.matrix <- rbind( 
  "Trisomies"=hrd.counts[, final.clean.samples],
  final.tx.numeric[tx.order[-1], final.clean.samples],
  final.snv.numeric.matrix[,final.clean.samples],
  final.cnv.numeric.matrix[,final.clean.samples], 
  final.sv.numeric.matrix[,final.clean.samples]
)

# Metadata annotation -----------------------------------------------------


## Set colors -------------------------------------------------------------


# add hrd as first row
n.tx <- nrow(final.tx.numeric)
n.snv <- nrow(final.snv.numeric.matrix)
n.cnv <- nrow(final.cnv.numeric.matrix); 
cnv.colors=ifelse( startsWith(rownames(final.cnv.numeric.matrix), "-")|startsWith(rownames(final.cnv.numeric.matrix), "del"), "#7570B3", "#E7298A")
n.sv <- nrow(final.sv.numeric.matrix)
row.names.colors <- c("#4DAF4A", rep("#377EB8", 3), "#E41A1C", rep("#984EA3", 3), "darkred", rep("#111111", n.snv), cnv.colors, rep("#B30000", n.sv))
row.names.fontface <- c(1, rep(1, n.tx), rep(3, n.snv), rep(1, length(cnv.colors)), rep(3, n.sv))
mutation.highlight <- c("HNRNPU", "FIP1L1", "IKZF3", "IKBKB", "HNRNPA2B1", "SGPP1", 
                        "SP3", 
                        # "+16p13.13", "+3q26.2", # Rustad et al has ≥ 100 peaks so pretty much impossible to "claim" new ones. Let's remove new.
                        "PRSS2 / MTRNR2L6", 
                        "PRR14L", "MFF", "EMSY", "SDCCAG8", "GPR180",
                        "ICE1")
to.highlight <- which(rownames(final.numeric.matrix) %in% mutation.highlight)
row.names.fontface[to.highlight] <- row.names.fontface[to.highlight]+1
imwg.colors <- c("Hyperdiploid"="#4DAF4A", "Cyclin D"="#377EB8", "MMSET"="#E41A1C", "MAF"="#984EA3", "Unclassified"="#999999")
stage.colors <- c("MGUS"= "#006D2C", "SMM"="#66C2A4", "MM"="#E5F5E0")
assay.colors <- c("WGS"= "#756BB1", "W(G+X)S"="#CBC9E2", "WXS"="#F2F0F7")
cohort.colors <- c("SU2C"="#386CB0", "CTC"="#F0027F", "Bustoros"="#BF5B17", "MMRF"="#EEEEEE", "OBEN"="#BEAED4", "WTC"="#FDC086")
gender.colors <- c("F"="#FCCDE5", "M"="#80B1D3")
age.colors <- circlize::colorRamp2(seq(min(meta.df$Age), max(meta.df$Age), length=9), brewer.pal(9, "BuPu"))

c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5")
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
  "#E5C494", "#B3B3B3")
c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
  "#A6761D", "#666666")

Gain.Palette <- rev(brewer.pal(10, "PiYG")[1:5])
# Loss.Palette <- brewer.pal(4, "Blues")
Loss.Palette <- brewer.pal(4, "Purples") #isolate pinks from green to pink gradient

class.values.colors <- c(
  "Negative"="white",
  
  "Gain1"= Gain.Palette[2],
  "Gain2"=Gain.Palette[3],
  "Gain3"=Gain.Palette[4],
  
  "Loss1"=Loss.Palette[2],
  "Loss2"=Loss.Palette[3], 
  "Loss3"=Loss.Palette[4],
  
  "HRD1"="white", 
  "HRD2"="#E5F5E0", 
  "HRD3"="#C7E9C0",
  "HRD4"="#A1D99B",
  "HRD5"="#74C476",
  "HRD6"="#4DAF4A",
  "HRD7"="#4DAF4A",
  "HRD8"="#4DAF4A",
  
  "SNV/INDEL"="#222222",
  "SV1"="#BEBADA",
  "SV2"="darkred",
  "t(11;14)(CCND1)"="#377EB8",
  "t(12;14)(CCND2)"="#377EB8",
  "t(14;16)(MAF)"="#984EA3",
  "t(14;20)(MAFB)"="#984EA3",
  "t(4;14)(MMSET)"="#E41A1C",
  "t(6;14)(CCND3)"="#377EB8",
  "t(8;14)(MAFA)"="#984EA3",
  "t(MYC;IgH/K/L)"="darkred")

patient.annotation.df <- meta.df[final.clean.samples, c("IMWG", "Stage", "Cohort", "Assay", "Age", "Gender")]

patient.annotation <- HeatmapAnnotation(df=patient.annotation.df[, c("IMWG", "Stage", "Age", "Gender")], 
                                        col = list(IMWG=imwg.colors, Stage=stage.colors, Cohort=cohort.colors, Age=age.colors, Assay=assay.colors, Gender=gender.colors),
                                        height = unit(1.5, "cm"),
                                        simple_anno_size_adjust = TRUE,
                                        annotation_name_gp= gpar(fontsize = 9))

# gene left annotation ----------------------------------------------------

gene.q.val.for.final.matrix <- gene.q.val[rownames(final.numeric.matrix),"q"]
gene.groups.for.final.matrix <- gene.q.val[rownames(final.numeric.matrix),"group"]

gene.q.ha <- rowAnnotation("Significance\n-10log10(FDR)"=anno_barplot(pmin(-10*log10(gene.q.val.for.final.matrix), 100), 
                                                                      axis_param = list(direction = "reverse")), 
                           annotation_name_gp = gpar(fontsize = 7))

# OR annotation -----------------------------------------------------------

# will do a simple left join, group by, nest a fisher test and hop!
all(colnames(final.numeric.matrix)==rownames(meta.df[final.clean.samples, ]))

# 1030 patients
# test.df contains patients for OR significant test
# For mutation, we will used powered genomes and validate in exomes
# For CNV and SV, we don't really have a model, so we'll use same as for mutations
# For initiating event, which is also annotated from patient charts, 
#    we will use all genomes regardless of power to detect the translocation 
#    (which is hard to estimate anyways)
test.df <- as.data.frame( cbind(select( meta.df[final.clean.samples,], "Stage"), t(final.numeric.matrix)) )

write_tsv(test.df |> tibble::rownames_to_column("Participant_ID"), "../data/fig1_all_matrix_patients_with_disease_stage.tsv")

all(colnames(test.df)[-1] == rownames(final.numeric.matrix))

# this to test difference in prevalence across MM vs MGUS for a given driver
all(clinical.and.terra.ref.pairs.samples$Participant_ID %in% rownames(test.df))

# remove exomes (will be used to validate signatures using significant genes - want to remove all sources of leakage)
exomes.smm.jco <- which(startsWith(rownames(test.df), "SMM"))
test.df <- test.df[-exomes.smm.jco, ]
dim(test.df)
# > dim(test.df)
# [1] 969 108
write_tsv(test.df|> tibble::rownames_to_column("Participant_ID"), "../data/fig1_all_matrix_patients_except_jco_with_disease_stage.tsv")


translocations <- rownames(gene.q.val[gene.q.val$group=="Translocation",])

# remove underpowered for testing odds ratio significant driver
underpowered.mgus.smm <- clinical.and.terra.ref.pairs.samples |> 
  filter(power_ccf_1 < 0.8) |> 
  pull(Participant_ID)

# mutations with weights
# say 1 for missense, 2 for nonsense; etc, 1 for subclonal CN, 2 for clonal, etc. 
# inspired by chapuy clustering
long.test.matrix.for.odds.ratio <- test.df %>%
  tibble::rownames_to_column("Sample ID") %>%
  pivot_longer(cols = -c(Stage, `Sample ID`))

# binarize for detection test
# keep underpowered for translocactions since forcecalling + FISH reports + initiating events

long.test.matrix.for.odds.ratio <- long.test.matrix.for.odds.ratio %>%
  mutate(value=case_when(name=="Trisomies" & value<2 ~ 0,
                         name=="Trisomies" & value>=2 ~ 1,
                         value>=1 ~ 1,
                         TRUE ~ 0)) %>%
  mutate(Stage=ifelse(Stage=="MM", "MM", "MGUS/SMM")) |>
  filter(name %in% translocations | ! ( `Sample ID` %in% underpowered.mgus.smm))


## Select powered samples --------------------------------------------------

## FISHER TESTS ------------------------------------------------------------
# to estimate enrichment in MM vs MGUS/SMM in well-powered individuals

fishers <- long.test.matrix.for.odds.ratio |> 
  nest(data = -name) %>%
  mutate(
    fit = map(data, ~ fisher.test(.x$value, .x$Stage)),
    tidied = map(fit, tidy)) %>%
  unnest(tidied) %>%
  select(-data, -fit)
fishers$group <- gene.q.val[fishers$name, "group"]

fishers <- fishers %>%
  group_by(group) %>%
  mutate(q.value=p.adjust(p.value, method="fdr")) %>%
  mutate(q.score=-10*log10(q.value))

write_tsv(fishers, "../data/fig1_fisher_tests_adjusted_pvalues.tsv", quote = "all")

## BINOM TESTS -------------------------------------------------------------
# to extract confidence intervals around observed frequencies in subgroups

binoms <- long.test.matrix.for.odds.ratio %>%
  nest(data = -name) %>%
  mutate(
    freq.mm.binom = map(data, ~binom_test(sum(.x$value==1 & .x$Stage=="MM", na.rm=TRUE), sum(.x$Stage=="MM" & !is.na(.x$value), na.rm=TRUE))),
    freq.smm.binom = map(data, ~binom_test(sum(.x$value==1 & .x$Stage=="MGUS/SMM", na.rm=TRUE), sum(.x$Stage=="MGUS/SMM" & !is.na(.x$value), na.rm=TRUE))),
  ) %>%
  unnest(freq.mm.binom, freq.smm.binom) %>%
  rename(estimate.mm=estimate, conf.low.mm=conf.low, conf.high.mm=conf.high, estimate.smm=estimate1, conf.low.smm=conf.low1, conf.high.smm=conf.high1) %>%
  select(name, estimate.mm, conf.low.mm, conf.high.mm, estimate.smm, conf.low.smm, conf.high.smm)

write_tsv(binoms, "../data/fig1_binom_confint_prevalence.tsv")

# supp table 4
fishers <- inner_join(fishers, binoms)

write_tsv(fishers, "../data/fig1_fisher_tests_adjusted_pvalues_with_binom.tsv", quote = "all")

# Row annotation ----------------------------------------------------------

or.df <- fishers %>% 
  tibble::column_to_rownames("name") %>% 
  mutate(estimate=ifelse(estimate%in%c(Inf, 0), 0, log(estimate))) %>%
  mutate(signif=ifelse(q.value<0.1, "*", " ")) %>%
  mutate(color=case_when(estimate>0 & q.value<0.1 ~ "#E78AC3", estimate<0 & q.value<0.1~"#006D2C", TRUE ~ "white")) %>%
  as.data.frame()
or.df <- or.df[rownames(final.numeric.matrix),]

or.col <- rowAnnotation("Association MM vs MGUS/SMM\n(log odds ratio)" = anno_barplot(select(or.df, "estimate"), 
                                                                                      baseline = 0, 
                                                                                      gp = gpar(fill = or.df[["color"]])
),
"S"=anno_text(or.df[,"signif"],  gp = gpar(fontsize = 7)),
annotation_name_gp = gpar(fontsize = 7))

# Actual heatmap ----------------------------------------------------------

set.seed(1234)

ht_opt(message=FALSE, 
       legend_labels_gp=gpar(fontsize=7),
       legend_title_gp=gpar(fontsize=7, fontface=2))

# get column order
Hm.numeric <- Heatmap(final.numeric.matrix, show_heatmap_legend = FALSE,
                      column_split=meta.df[final.clean.samples, "IMWG"],
                      column_title_gp = gpar(fontsize = 7),
                      row_split=gene.q.val[rownames(final.numeric.matrix), "group"],
                      row_title_gp = gpar(fontsize = 7),
                      show_column_names = FALSE,
                      cluster_rows = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      col = c("white", "black"),
                      row_names_gp = gpar(fontsize = 5, col=row.names.colors, fontface=row.names.fontface),
                      border = TRUE)
hm_col_order <- column_order(Hm.numeric)


# Class heatmap instead of numeric ----------------------------------------

# after discussions - this to have more freedom in coloring events int he actual figure
# hierarchical clustering is not used anyway; instead we prefered a classical heatmap

final.class.matrix <- final.numeric.matrix

# convert numeric values to class
table(gene.q.val[rownames(final.numeric.matrix), "group"])

# Translocation
table(rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Translocation"])
table(final.class.matrix[gene.q.val[rownames(final.numeric.matrix), "group"]=="Translocation"])

final.class.matrix["Trisomies",] <- paste0("HRD", final.class.matrix["Trisomies",])
other.translocations <- c("t(11;14)(CCND1)", "t(6;14)(CCND3)", "t(12;14)(CCND2)", "t(4;14)(MMSET)", 
                          "t(14;16)(MAF)", "t(14;20)(MAFB)", "t(8;14)(MAFA)", "t(MYC;IgH/K/L)")
final.class.matrix[other.translocations,] <- ifelse(final.class.matrix[other.translocations,] >= 1, rownames(final.class.matrix[other.translocations,]), "Negative")
final.class.matrix[final.class.matrix%in%c("HRD0")] <- "Negative"

# Point mutation\nSmall indel
snvs <- rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Point mutation\nSmall indel"]
table(rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Point mutation\nSmall indel"])
table(final.class.matrix[gene.q.val[rownames(final.numeric.matrix), "group"]=="Point mutation\nSmall indel"])
final.class.matrix[snvs,] <- ifelse(final.class.matrix[snvs,] >= 1, "SNV/INDEL", "Negative")

# Copy number\nAbnormality
cnvs <- rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Copy number\nAbnormality"]
table(rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Copy number\nAbnormality"])
table(final.class.matrix[gene.q.val[rownames(final.numeric.matrix), "group"]=="Copy number\nAbnormality"])
cnv.gains <- cnvs[startsWith(cnvs, "+")]
cnv.losses <- cnvs[startsWith(cnvs, "-")]
final.class.matrix[cnv.gains,] <- paste0("Gain", final.class.matrix[cnv.gains,])
final.class.matrix[cnv.losses,] <- paste0("Loss", final.class.matrix[cnv.losses,])
final.class.matrix[final.class.matrix%in%c("Gain0", "Loss0")] <- "Negative"

# Structural\nVariant
table(final.numeric.matrix[gene.q.val[rownames(final.numeric.matrix), "group"]=="Structural\nVariant"])
svs <- rownames(final.numeric.matrix)[gene.q.val[rownames(final.numeric.matrix), "group"]=="Structural\nVariant"]
final.class.matrix[svs,] <- paste0("SV", final.class.matrix[svs,]) # 1 no LOF, 2 LOF
final.class.matrix[final.class.matrix%in%c("SV0")] <- "Negative"

table(final.class.matrix, useNA = "ifany")

place.n.samples = rep("", ncol(final.class.matrix))
place.n.samples[round(length(place.n.samples)/2)] <- paste0("N = ", scales::comma(ncol(final.class.matrix)))
n.samples = HeatmapAnnotation(n_samples = anno_text(place.n.samples, location = 1, rot = 0, just = "right", gp = gpar(fontsize = 7)))


write.table(final.class.matrix, "../data/figure1_matrix.tsv", sep="\t")
write.table(patient.annotation.df, "../data/figure1_annot.tsv", sep = "\t")


# Summary stats -----------------------------------------------------------

dim(patient.annotation.df)
table(patient.annotation.df$Stage)

# how many samples
patient.annotation.df |> 
  filter(Cohort %in% c("CTC", "SU2C")) |>
  group_by(Stage) |>
  count()

# how many samples
patient.annotation.df |> 
  filter(Cohort %in% c("CTC", "SU2C")) |>
  group_by(Stage, Cohort) |>
  count()

patient.annotation.df |> 
  filter(Cohort %in% c("MMRF")) |>
  group_by(Stage, Cohort) |>
  count()

# how many MGIS/SMM genomes total
patient.annotation.df |>
  filter(Assay=="WGS" & Stage %in% c("MGUS", "SMM")) |>
  group_by(Stage) |>
  count()

# how many have a driver
table( colSums( final.class.matrix[snvs, patient.annotation.df |> filter(Cohort %in% c("CTC", "SU2C") & Stage %in% c("MGUS", "SMM")) |> rownames()] != "Negative" ))

# how many drivers in MGUS
summary( colSums( final.class.matrix[c(snvs), patient.annotation.df |> filter(Stage %in% c("MGUS")) |> rownames()] != "Negative" ))
summary( colSums( final.class.matrix[c(snvs), patient.annotation.df |> filter(Stage %in% c("SMM")) |> rownames()] != "Negative" ))
summary( colSums( final.class.matrix[c(snvs), patient.annotation.df |> filter(Stage %in% c("MM")) |> rownames()] != "Negative" ))

mgus.smm.no.driver <- colSums( final.class.matrix[c(snvs), patient.annotation.df |> rownames()] != "Negative" ) == 0

sort( colSums( final.class.matrix[c(svs, cnvs, tx.order), mgus.smm.no.driver] != "Negative"))
sort( colSums( final.class.matrix[c(tx.order), mgus.smm.no.driver] != "Negative"))

summary( colSums( final.class.matrix[c(cnvs), patient.annotation.df |> filter(Stage %in% c("MM")) |> rownames()] != "Negative" ))

# how many have SV drivers or IgH translocation
summary( colSums( final.class.matrix[c(tx.order[-1], svs), patient.annotation.df |> filter(Stage %in% c("MM")) |> rownames()] != "Negative" ))
summary( colSums( final.class.matrix[c(tx.order[-1], svs), patient.annotation.df |> filter(Stage %in% c("MGUS")) |> rownames()] != "Negative" ))
# the effect of SV hotspots
table( colSums( final.class.matrix[c(tx.order[-1], svs), patient.annotation.df |> filter(Assay != c("WXS")) |> rownames()] != "Negative" ) >= 1) / patient.annotation.df |> filter(Assay != c("WXS")) |> nrow()
# IgH hotspots
table( colSums( final.class.matrix[c(tx.order[2:8], svs), patient.annotation.df |> filter(Assay != c("WXS")) |> rownames()] != "Negative" ) >= 1) / patient.annotation.df |> filter(Assay != c("WXS")) |> nrow()
# No IgH CatchTheFISH But SV variants
table(colSums(final.class.matrix[tx.order[2:7],]!="Negative")==0 & final.class.matrix["NSD2",]!="Negative") # 9 to investigate
table(colSums(final.class.matrix[tx.order[2:7],]!="Negative")==0 & final.class.matrix["WWOX (MAF)",]!="Negative") # 14 to investigate

non.igh.svs.hotspots <- svs[-which(svs%in%c("WWOX (MAF)", "NSD2"))]

table(colSums(final.class.matrix[non.igh.svs.hotspots, patient.annotation.df |> filter(Assay != c("WXS") & Stage == "MM") |> rownames()] != "Negative"))
prop.test(218, 812) # MM with a SV, p=0.2684729, CI 0.2385431 0.3006270
table(colSums(final.class.matrix[non.igh.svs.hotspots, patient.annotation.df |> filter(Assay != c("WXS") & Stage != "MM") |> rownames()] != "Negative"))
prop.test(12, 157) # MGUS/SMM with a SV, p=0.07643312, CI  0.04190557 0.13268084
fisher.test(matrix(c(218, 812, 12, 157), nrow=2))
# p-value = 4.323e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.906809 7.070659
# sample estimates:
#   odds ratio 
# 3.509516 

# > non.igh.svs.hotspots
# [1] "TRAF3 / AMN"      "SP140 / SP140L"   "PRSS2 / MTRNR2L6" "CDKN2C"           "PRR14L"           "MFF"              "EMSY"             "RB1 SV"           "SDCCAG8"         
# [10] "CDKN2A"           "GPR180"           "ICE1"             "TBL1XR1" 


# HEATMAP CALL ------------------------------------------------------------

# March 25 split by initiating event, then reorder by stage and then mutations
# > pivot_longer(-`Sample ID`, names_to = "Mutation_Name", values_to = "Mutation_Status")

do.not.order.columns.based.on <- c("Trisomies", "t(MYC;IgH/K/L)")
all.drivers.to.order.on <- rownames(gene.q.val[rownames(final.numeric.matrix)[-which(rownames(final.numeric.matrix) %in% do.not.order.columns.based.on)],])

variable.object.to.order.on <- c("Stage", all.drivers.to.order.on)
numeric.matrix.to.order.samples <- (-1*(final.numeric.matrix>0)) |> t() |> as.data.frame() |> rownames_to_column("Sample ID")

new.order <- meta.df[colnames(final.numeric.matrix),] |> 
  rownames_to_column("Sample ID") |>
  left_join(numeric.matrix.to.order.samples) |>
  arrange(!!! rlang::syms(variable.object.to.order.on)) |>
  pull(`Sample ID`)
# > new.order[new.order %nin% colnames(final.class.matrix)]
# [1] "MMRF_1889" "MMRF_2363" # these two samples are excluded (not the baseline sample if I recall correctly)

Hm <- Heatmap(final.class.matrix, 
              column_order = new.order,
              show_heatmap_legend = FALSE,
              column_split=meta.df[final.clean.samples, c("IMWG")],
              column_title_gp = gpar(fontsize = 7),
              row_split=gene.q.val[rownames(final.numeric.matrix), "group"],
              row_title_gp = gpar(fontsize = 7),
              show_column_names = FALSE,
              cluster_rows = FALSE,
              cluster_column_slices = FALSE,
              show_column_dend = FALSE,
              col = class.values.colors,
              row_names_gp = gpar(fontsize = 5, col=row.names.colors, fontface=row.names.fontface),
              border = TRUE,
              left_annotation = gene.q.ha,
              right_annotation = or.col, 
              bottom_annotation = n.samples)

ht.list <- patient.annotation %v% Hm

draw(ht.list)

pdf("../figures/FIG1_Mutations_with_class_stage.pdf", width = 6, height = 8)
draw(ht.list)
dev.off()

