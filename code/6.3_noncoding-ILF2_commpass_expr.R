setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(DESeq2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(RColorBrewer)

expr.commpass <- readRDS("../Andrea/SU2C_CoMM/results/vst_normalized_gene_expression_counts_CoMMpas.rds")
class(expr.commpass)
dim(expr.commpass)
colnames(expr.commpass)
rownames(expr.commpass)

# (This is the final list from extended data figure 12)
pts.ilf2.final <- c("MMRF_1290_1_BM_CD138pos", "MMRF_1389_1_BM_CD138pos", "MMRF_1403_1_BM_CD138pos", 
              "MMRF_1491_1_BM_CD138pos", "MMRF_1683_1_BM_CD138pos", "MMRF_1796_1_BM_CD138pos", 
              "MMRF_1855_1_BM_CD138pos", "MMRF_1912_1_PB_CD138pos", "MMRF_1965_1_PB_CD138pos", 
              "MMRF_2076_1_BM_CD138pos", "MMRF_2091_1_BM_CD138pos", "MMRF_2200_1_BM_CD138pos", 
              "MMRF_2272_1_BM_CD138pos", "MMRF_2535_1_BM_CD138pos", "MMRF_2739_1_BM_CD138pos", 
              "MMRF_2776_1_BM_CD138pos")
length(pts.ilf2.final)
# 16

# Mut ---------------------------------------------------------------------

#load pileups from small wolF tool
mutant <- read_tsv("../data/MPILEUP_153671247.txt", col_names = FALSE)
patients <- read_tsv("../data/PATIENTS", col_names = "PATIENT")
patients$pair <- basename(patients$PATIENT)
isSorted(patients$pair)

colnames.mutant <- c("Chromosome", "Position", "Reference", 
                     sort(c(paste0(patients$pair, "_DEPTH"),
                         paste0(patients$pair, "_SEQ"),
                         paste0(patients$pair, "_SEQ_QUAL"))))
colnames(mutant) <- colnames.mutant

MUT_COUNT=mutant %>% 
  mutate_all(as.character) %>%
  pivot_longer(cols=-c(1:3), 
               names_to = c("PAIR", "x"), 
               names_pattern = "(.*)_(DEPTH|SEQ|SEQ_QUAL)$") %>%
  pivot_wider(id_cols = 1:4, names_from = x, values_from = value) %>%
  mutate(MUTANT_SEQ=str_remove_all( toupper(SEQ), "G"), N_MUTANT=nchar(MUTANT_SEQ))


# Forcecall ILF2 mutants -------------------------------------------

length(unique(MUT_COUNT$PAIR))
MUT_COUNT$PARTICIPANT_ID=str_extract(MUT_COUNT$PAIR, "MMRF_[0-9]{4,}")

length(unique(MUT_COUNT$PAIR))
length(unique(MUT_COUNT$PARTICIPANT_ID))

mm <- read_tsv("../data/export_MM-like_features_MGUS_SMM_MM.tsv")
MUT_COUNT <- MUT_COUNT |> inner_join(mm, by=c("PARTICIPANT_ID"="Sample ID"))
table(unique(MUT_COUNT$PARTICIPANT_ID) %in% mm$`Sample ID`)

mut.counted <- MUT_COUNT |> 
  mutate(Tumor_Seq_Allele2=str_extract(MUTANT_SEQ, "[ACT]"), 
         t_alt_count=str_count(MUTANT_SEQ, Tumor_Seq_Allele2),
         t_ref_count=str_count(toupper(SEQ), "G")) |>
  mutate(Tumor_Seq_Allele2=replace_na(Tumor_Seq_Allele2, "G"),
         t_alt_count=replace_na(t_alt_count, 0)) |>
  mutate(expr.pair= str_remove(PAIR, "_pair"))

mut.counted |>
  filter((t_alt_count+t_ref_count)>=2) |>
  select(PARTICIPANT_ID) |>
  distinct() |>
  nrow()
# 770 participants from CoMMpass with at least 2 reads on the hotspot mutations (either or)

mut.counted.hotspot <- mut.counted |>
  mutate(Pos=paste0(Chromosome, Position)) |>
  select(PAIR, PARTICIPANT_ID, expr.pair, amp_chr_1q, amp_G.peak.1_amp_1q21.2, Gender, Age, Stage, HRDTx, Cohort, IMWG, Chromosome, Pos, t_ref_count, t_alt_count) |>
  pivot_wider(id_cols = c(PAIR, PARTICIPANT_ID, Gender, Age, Stage, HRDTx, Cohort, IMWG, Chromosome, expr.pair, amp_chr_1q, amp_G.peak.1_amp_1q21.2,), 
              names_from = Pos, 
              values_from = c(t_ref_count, t_alt_count)) |>
  arrange(PAIR) |>
  group_by(PARTICIPANT_ID) |>
  slice_head(n=1)

length(unique(mut.counted.hotspot$PARTICIPANT_ID))
length(unique(mut.counted.hotspot$PAIR))


# Classify well-covered and mutants ---------------------------------------

mut.counted.hotspot <- mut.counted.hotspot |>
  mutate(COVERED = ((t_ref_count_chr1153671247+t_alt_count_chr1153671247>=2)|(t_ref_count_chr1153671254+t_alt_count_chr1153671254>=2))) |>
  mutate(MUTANT = (COVERED & ( t_alt_count_chr1153671247>=1 | t_alt_count_chr1153671254>=1))) |>
  mutate(N_MUT = case_when(MUTANT==FALSE ~ 0, 
                           MUTANT==TRUE &  t_alt_count_chr1153671247>=1  & t_alt_count_chr1153671254>=1 ~ 2,
                           MUTANT==TRUE ~ 1,
                           TRUE ~ -1)) |>
  mutate(ILF2_MUT=MUTANT,
         Gain1q= amp_chr_1q >=1 | amp_G.peak.1_amp_1q21.2 >=1,
         group = case_when( 
           ILF2_MUT & Gain1q ~ "ILF2 MUT\n+1q21.2",
           ILF2_MUT & !Gain1q ~ "ILF2 MUT",
           Gain1q & !ILF2_MUT ~ "ILF2 WT\n+1q21.2",
           !Gain1q & !ILF2_MUT ~ "ILF2 WT",
           TRUE ~ "WHAT"))

table(mut.counted.hotspot$COVERED)
# FALSE  TRUE 
# 23   769 
table(mut.counted.hotspot$MUTANT)
# FALSE  TRUE 
# 766   26 
table(mut.counted.hotspot$MUTANT, mut.counted.hotspot$IMWG)
fisher.test(table(mut.counted.hotspot$MUTANT, mut.counted.hotspot$IMWG=="MAF"))
# Fisher's Exact Test for Count Data
# 
# data:  table(mut.counted.hotspot$MUTANT, mut.counted.hotspot$IMWG == "MAF")
# p-value = 4.702e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   9.176794 60.135566
# sample estimates:
# odds ratio 
#    23.4112 

table(mut.counted.hotspot$MUTANT, mut.counted.hotspot$N_MUT)

colnames(mut.counted.hotspot)

table.s12 <- mut.counted.hotspot |>
  filter(MUTANT==TRUE) |>
  select(PARTICIPANT_ID, PAIR, IMWG, amp_chr_1q, amp_G.peak.1_amp_1q21.2,
         "t_ref_count_chr1153671247", "t_alt_count_chr1153671247", "t_ref_count_chr1153671254","t_alt_count_chr1153671254", 
         N_MUT)
nrow(table.s12)
table.s12 |> write_tsv("../data/tableS12_ILF2mutants.tsv")
table.s12 |> clipr::write_clip()

mut.counted.hotspot |>
  filter(MUTANT==TRUE) |>
  pull(PARTICIPANT_ID)
# [1] "MMRF_1328" "MMRF_1389" "MMRF_1403" "MMRF_1491" "MMRF_1602" "MMRF_1618" "MMRF_1683" "MMRF_1796" "MMRF_1816" "MMRF_1855"
# [11] "MMRF_1986" "MMRF_2041" "MMRF_2076" "MMRF_2091" "MMRF_2153" "MMRF_2185" "MMRF_2200" "MMRF_2226" "MMRF_2272" "MMRF_2331"
# [21] "MMRF_2535" "MMRF_2637" "MMRF_2739" "MMRF_2743" "MMRF_2761" "MMRF_2776"
large.ilf2.mutants <- mut.counted.hotspot |>
  filter(MUTANT==TRUE) |>
  pull(PARTICIPANT_ID)

mut.counted.hotspot

table(mut.counted$expr.pair %in% colnames(expr.commpass))
table(colnames(expr.commpass) %in% mut.counted$expr.pair)

colnames(expr.commpass)[!(colnames(expr.commpass) %in% mut.counted$expr.pair)]

# 1Mb around
# also GM12878 TAD
# chr1:153120000-154200000


# Extract RNA expression in neighborhood ------------------------------------------------------------

window.start = 152652524
window.end = 154474524

expr.df <- expr.commpass %>% 
  # filter(hgnc_symbol=="ILF2") %>%
  filter(nchar(hgnc_symbol)>1 & chromosome_name=="1" & start_position >= 	window.start & end_position <= window.end) |>
  pivot_longer(starts_with("MMRF_"), names_to = "expr.pair", values_to = "vst_norm")

# export list of genes
expr.df |> 
  select(1:6) |> 
  distinct() |>
  write_tsv("../data/ilf2-neighbor-tads.genes.tsv")


# Join Expression and Mutation status -------------------------------------------------------

mut.expr.df <- expr.df |>
  inner_join(mut.counted, by=c("expr.pair"))

mut.expr.df.any.mutant <- expr.df |>
  inner_join(mut.counted.hotspot, by=c("expr.pair"))

# how many covered and expr
mut.expr.df.any.mutant.covered <- mut.expr.df.any.mutant |> 
  filter(COVERED==TRUE)
mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2") |> group_by(COVERED) |> count()
nrow(mut.expr.df.any.mutant.covered|> filter(hgnc_symbol=="ILF2"))

# Annotate Group ---------------------------------------------------------------
# for single position only; any hotspot is already merged with matrix
meta.final <- read_tsv("../data/FINAL_meta_columns_matrix.tsv")

mut.expr.df <- mut.expr.df |>
  mutate(MMRF=str_extract(expr.pair, "MMRF_[0-9]{4}")) |>
  left_join(meta.final, by=c("MMRF"="Sample ID"))

cn.final <- read_tsv("../data/FINAL_Nov2_CNV_Focal_Numeric_Matrix.tsv")

cn.final <- cn.final |>
  filter(Variable%in%c("+1q", "+1q21.2")) |>
  pivot_longer(-Variable, names_to = "Sample ID", values_to = "CN") |>
  pivot_wider(id_cols = "Sample ID", names_from = "Variable", values_from = "CN")

table(cn.final$`Sample ID` %in% mut.expr.df$MMRF)
table(mut.expr.df$MMRF %in% cn.final$`Sample ID`)

mut.expr.df <- mut.expr.df |>
  left_join(cn.final, by=c("MMRF"="Sample ID"))

mut.expr.df <- mut.expr.df |> 
  
  filter( (t_alt_count+t_ref_count) >= 2 ) |>
  
  mutate(ILF2_MUT=t_alt_count>=1) |>
  
  mutate(group = case_when( 
    t_alt_count>=1 & `+1q21.2` >=1 ~ "ILF2 MUT\n+1q21.2",
    t_alt_count>=1 ~ "ILF2 MUT",
    `+1q21.2` >=1 ~ "ILF2 WT\n+1q21.2",
    TRUE ~ "ILF2 WT"
                            )) |>
  mutate(group = fct_reorder(group, vst_norm, median) )

mut.expr.df |> 
  filter(hgnc_symbol=="ILF2") |>
  group_by(group) |>
  count()

mut.expr.df.any.mutant.covered <- mut.expr.df.any.mutant.covered |>
  filter(COVERED==TRUE) |>
  mutate(ILF2_MUT=MUTANT,
         Gain1q= amp_chr_1q >=1 | amp_G.peak.1_amp_1q21.2 >=1,
         group = case_when( 
           ILF2_MUT & Gain1q ~ "ILF2 MUT\n+1q21.2",
           ILF2_MUT & !Gain1q ~ "ILF2 MUT",
           Gain1q & !ILF2_MUT ~ "ILF2 WT\n+1q21.2",
           !Gain1q & !ILF2_MUT ~ "ILF2 WT",
           TRUE ~ "WHAT"))

nrow(mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2"))
# > mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2") |> group_by(group) |> count()
# # A tibble: 4 × 2
# # Groups:   group [4]
# group                   n
# <chr>               <int>
#   1 "ILF2 MUT"             12
# 2 "ILF2 MUT\n+1q21.2"     9
# 3 "ILF2 WT"             330
# 4 "ILF2 WT\n+1q21.2"    238
mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2") |> nrow()

# Stats -------------------------------------------------------------------

# 
stats.any.ILF2 <- mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2") |> do(tidy(lm(vst_norm~ILF2_MUT+Gain1q, data=.)))
stats.any.ILF2
# # A tibble: 3 × 5
# term         estimate std.error statistic  p.value
# <chr>           <dbl>     <dbl>     <dbl>    <dbl>
#   1 (Intercept)    12.3      0.0241    511.   0       
# 2 ILF2_MUTTRUE    0.233    0.0981      2.38 1.78e- 2
# 3 Gain1qTRUE      0.551    0.0369     15.0  4.33e-43

stats.any.ILF2.single.double <- mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="ILF2") |> do(tidy(lm(vst_norm~N_MUT+Gain1q, data=.)))
stats.any.ILF2.single.double

stats.any.TPM3 <- mut.expr.df.any.mutant.covered |> filter(hgnc_symbol=="TPM3") |> do(tidy(lm(vst_norm~ILF2_MUT+Gain1q, data=.)))
stats.any.TPM3

stats.any <- mut.expr.df.any.mutant.covered |> group_by(hgnc_symbol) |> do(tidy(lm(vst_norm~ILF2_MUT+Gain1q, data=.)))
stats.any <- stats.any |> group_by(term) |> mutate(q=p.adjust(p.value, method="fdr")) |> ungroup()
stats.any.ilf2.term <- stats.any |> filter(term=="ILF2_MUTTRUE") |> arrange(p.value)
# 
# ### following code splits per mutation - not used (each of the two nucleotide)
# fig mut 1
# counts.pos.1 <- mut.expr.df |>
#   filter(Position==153671247) |>
#   group_by(group, hgnc_symbol) |>
#   summarise(n=n(), label=paste0("italic('N=", n, "')"))
# 
# #  (hg38) 153671247 
# #  (hg19) 153643723
# 
# expr.ilf2.1 <- ggplot(mut.expr.df |> filter(Position==153671247 & hgnc_symbol=="ILF2"), aes(group, vst_norm, color=group)) +
#   geom_boxplot() +
#   geom_text(data=counts.pos.1, aes(label=label), y=11, color="black", size=1.5, parse = TRUE) +
#   # geom_jitter(alpha=0.2) +
#   # stat_compare_means(method="anova", label="p", label.x.npc = "right", label.y.npc = "bottom", size=1.5) +
#   # stat_pvalue_manual(stats.pos.1.ILF2, label = "p.adj.signif", hjust = 1, hide.ns = TRUE, tip.length = 0, step.increase = 0.05) +
#   scale_color_manual(values=brewer.pal(11, "PiYG")[c(11, 3, 2, 1)]) +
#   # scale_y_continuous(limits=c(11, 15.5)) +
#   labs(x="", y="ILF2 Gene Expression\n(VST normalized value)", title = "chr1:153,643,723") +
#   theme_bw(base_size = 6, base_line_size = 0.25, base_rect_size = 0.25) +
#   theme(panel.grid = element_blank(), legend.position = "none")
#   # facet_wrap(~hgnc_symbol, scales = "free_y")
# 
# expr.ilf2.1
# 
# ggsave("../figures/commpass_ilf2.1.pdf", plot=expr.ilf2.1, width = 2, height = 1.5)
# 
# fig mut 2
# 
# counts.pos.2 <- mut.expr.df |>
#   filter(Position==153671254) |>
#   group_by(group, hgnc_symbol) |>
#   summarise(n=n(), label=paste0("italic('N=", n, "')"))
# 
# #  (hg38) 153671254
# #  (hg19) 153643730
# 
# expr.ilf2.2 <- ggplot(mut.expr.df |> 
#                         filter(Position==153671254 & hgnc_symbol=="ILF2"), aes(group, vst_norm, color=group)) +
#   geom_boxplot() +
#   geom_text(data=counts.pos.2, aes(label=label), y=11, color="black", size=1.5, parse = TRUE) +
#   # geom_jitter(alpha=0.2) +
#   # stat_compare_means(method="anova", label="p", label.x.npc = "right", label.y.npc = "bottom", size=1.5) +
#   # stat_pvalue_manual(stats.pos.2, label = "p.adj.signif", hjust = 1, hide.ns = TRUE, tip.length = 0, step.increase = 0.05) +
#   scale_color_manual(values=brewer.pal(11, "PiYG")[c(11, 3, 2, 1)]) +
#   scale_y_continuous(limits=c(11, 15.5)) +
#   labs(x="", y="ILF2 Gene Expression\n(VST normalized value)", title = "chr1:153,643,730") +
#   theme_bw(base_size = 6, base_line_size = 0.25, base_rect_size = 0.25) +
#   theme(panel.grid = element_blank(), legend.position = "none")
# 
# expr.ilf2.2
# 
# ggsave("../figures/commpass_ilf2.2.pdf", plot=expr.ilf2.2, width = 2, height = 1.5)


# export raw data

mut.expr.df |>
  select(hgnc_symbol, ensembl_gene_id, expr.pair, Chromosome, Position, Tumor_Seq_Allele2, t_alt_count, t_ref_count, vst_norm, `+1q`, `+1q21.2`) |>
  write_tsv("../data/ILF2_commpass_mut_counts.tsv")


# Neighbor gene for any mutant --------------------------------------------

# simple model
# vst_norm ~ ILF2_MUT+Gain1q


lm.values <- mut.expr.df.any.mutant.covered |>
  group_by(hgnc_symbol) |>
  do(tidy(lm(vst_norm ~ ILF2_MUT+Gain1q, .))) |>
  filter(term=="ILF2_MUTTRUE") |>
  ungroup() |>
  mutate(q=p.adjust(p.value, method="fdr"))

# from S100A6 to RAB13
# on the figure we show more than the tad but we also want to adjust p just within the TAD
lm.values.within.TAD <- mut.expr.df.any.mutant.covered |>
  filter(start_position>=153534599 & end_position<=153986358) |>
  group_by(hgnc_symbol) |>
  do(tidy(lm(vst_norm ~ ILF2_MUT+Gain1q, .))) |>
  filter(term=="ILF2_MUTTRUE") |>
  ungroup() |>
  mutate(q=p.adjust(p.value, method="fdr")) |>
  mutate(q=case_when(q<0.1~q, TRUE ~NA)) |>
  mutate(text_only_signif=case_when(q<0.1~hgnc_symbol, TRUE ~"")) |>
  mutate(text_not_signif=case_when(q<0.1~"", TRUE ~hgnc_symbol))

lm.values.within.TAD |> write_tsv("../data/GM12878_TAD_ILF2_neigbors_expr.tsv")
# lm.values.within.TAD |> filter(q<0.1) |> clipr::write_clip()

lm.values.with.coordinates <- lm.values |>
  left_join(expr.commpass)

lm.values.with.coordinates.ILF2 <- lm.values.with.coordinates |> 

  mutate(q=case_when(q<0.1~q, TRUE ~NA)) |>
  mutate(text_only_signif=case_when(q<0.02~hgnc_symbol, TRUE ~"")) |>
  mutate(text_not_signif=case_when(q<0.02~"", TRUE ~hgnc_symbol))


ilf2.correlations.any.hotspot <- ggplot(lm.values.with.coordinates.ILF2) +
  
  geom_vline(xintercept=c(153671247, 153671254), color="grey", linetype=2) +
  annotate(geom = "text", x=153848000, y=-0.7, label="ILF2 Promoter" ,color="lightgrey") +
  
  geom_rect(aes(xmin=start_position, xmax=end_position, ymin=0, ymax=estimate, fill=-log10(q)), alpha=1, color="black", linewidth=0.25) +
  geom_text(aes(x=1/2*(start_position+end_position), y=estimate, label=text_only_signif), nudge_x = 5e4, nudge_y = 0.05, color="royalblue") +
  # geom_text(aes(x=1/2*(start_position+end_position), y=estimate, label=text_not_signif), nudge_x = 5e4, nudge_y = 0.05, size=2.5) +
  
  
  scale_fill_distiller(na.value = "white", palette = "YlGnBu") +
  scale_x_continuous(labels = scales::comma) +
  
  labs(x="Genomic coordinate (chromosome 1)", y="ILF2 term estimate \nExpression ~ ILF2 + 1q status") +
  
  theme_classic2()

ilf2.correlations.any.hotspot

# 800x500pts
ggsave(filename = "../figures/ilf2_correlations_tads_any_hotspot.pdf", 
       plot = ilf2.correlations.any.hotspot, 
       width=7, height = 4, units = "in")

# manhattan-plot for correlations with significance as y axis ---------------------

ilf2.manhattan.any.hotspot <- ggplot(lm.values.with.coordinates.ILF2) +
  
  geom_vline(xintercept=c(153671247, 153671254), color="grey", linetype=2) +
  geom_hline(yintercept = -log10(0.05), color="grey") +
  
  annotate(geom = "text", x=153848000, y=6, label="ILF2 Promoter" ,color="lightgrey") +
  annotate(geom = "text", x=152848000, y=1.5, label="-log10(0.05)" ,color="lightgrey") +
  
  geom_point(aes(x=(end_position+start_position)/2, y=-log10(p.value), fill=estimate), alpha=1, shape=21, size=2) +
  
  geom_segment(aes(x=(end_position+start_position)/2,xend=(end_position+start_position)/2, 
                   y=0, yend=-log10(p.value)-0.1), alpha=1, color="black", linewidth=0.25) +
  geom_text(aes(x=1/2*(start_position+end_position), y=-log10(p.value), label=text_only_signif), nudge_x = 1E5, nudge_y = 0.25, color="royalblue") +
  # geom_text(aes(x=1/2*(start_position+end_position), y=estimate, label=text_not_signif), nudge_x = 5e4, nudge_y = 0.05, size=2.5) +
  
  scale_fill_distiller(na.value = "white", palette = "RdBu", ) +
  scale_color_distiller(na.value = "white", palette = "RdBu") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0, 0.1))) +
  
  labs(x="Genomic coordinate (chromosome 1)", 
       y="-log10(p) for ILF2 term\nGene Expression ~ ILF2 mutant + 1q status",
       fill="Beta term\n(ILF2 mut.)",
       color="Beta term\n(ILF2 mut.)") +
  
  theme_classic2()

ilf2.manhattan.any.hotspot

# 800x500pts
ggsave(filename = "../figures/ilf2.manhattan.any.hotspot.pdf", 
       plot = ilf2.manhattan.any.hotspot, 
       width=7, height = 4, units = "in")


# Volcano TAD -------------------------------------------------------------

volcano.within.tad <- ggplot(lm.values.within.TAD, aes(estimate, -log10(p.value))) +
  
  geom_vline(xintercept = 0, color="lightgrey") +
  geom_point(aes(color=q<0.1)) +
  geom_text(aes(label=text_only_signif, hjust=1, fontface=3), position = position_nudge(x = -0.06), size=2.5) +
  
  scale_x_continuous(limits=c(-1.2, 1.2)) +
  
  labs(x="Beta estimate for ILF2 mutant\nExpression ~ ILF2 Promoter mutant + Gain(1q)",
       y="-log10(p)") +
  
  theme_bw(base_size = 6, base_line_size = 0.25, base_rect_size = 0.25) +
  theme(panel.grid = element_blank(),
        legend.position = "none")

volcano.within.tad
ggsave("../figures/commpass_volcano_within_tad_ilf2_mutant.pdf", volcano.within.tad, width = 2, height = 1.5)


# ILF2 higher expr Boxplot --------------------------------------------------------

ILF2.ILF2mut.Gain1q.Regression <- mut.expr.df.any.mutant.covered |>
  filter(hgnc_symbol=="ILF2") |>
  do(tidy(lm(vst_norm ~ ILF2_MUT+Gain1q, .)))

mut.expr.df.any.mutant.covered |>
  filter(hgnc_symbol=="ILF2") |>
  mutate(group=factor(group, levels=c("ILF2 WT", "ILF2 MUT", "ILF2 WT\n+1q21.2", "ILF2 MUT\n+1q21.2"))) |>
  do(tidy(lm(vst_norm ~ group, .)))

# need pairwise comparison for boxplots
pairwise.comparisons <- mut.expr.df.any.mutant.covered |>
  filter(hgnc_symbol=="ILF2") |>
  mutate(group=factor(group, levels=c("ILF2 WT", "ILF2 MUT", "ILF2 WT\n+1q21.2", "ILF2 MUT\n+1q21.2"))) |>
  pairwise_t_test(vst_norm ~ group, 
                  p.adjust.method = "fdr", 
                  alternative="greater", 
                  comparisons = list(c("ILF2 WT", "ILF2 MUT"), 
                                     c("ILF2 WT", "ILF2 WT\n+1q21.2"), 
                                     c("ILF2 MUT", "ILF2 WT\n+1q21.2"),
                                     c("ILF2 WT\n+1q21.2", "ILF2 MUT\n+1q21.2")))

group.counts.ilf2.group <- mut.expr.df.any.mutant.covered |>
  filter(hgnc_symbol=="ILF2") |>
  mutate(group=factor(group, levels=c("ILF2 WT", "ILF2 MUT", "ILF2 WT\n+1q21.2", "ILF2 MUT\n+1q21.2"))) |>
  group_by(group) |>
  count() |>
  mutate(label=paste0("N=", n))

boxplot.t.test.pairwise.ILF2.Gain1q <- mut.expr.df.any.mutant.covered |>
  filter(hgnc_symbol=="ILF2") |>
  mutate(group=factor(group, levels=c("ILF2 WT", "ILF2 MUT", "ILF2 WT\n+1q21.2", "ILF2 MUT\n+1q21.2"))) |>
  ggplot(aes(group, vst_norm, color=group)) +
  
  geom_boxplot(linewidth=0.25, outlier.size = 0.25) +
  
  geom_jitter(alpha=0.2, size=0.25, shape=1) +
  
  geom_text(data=group.counts.ilf2.group, aes(y=11, label=label), color="#333333", size=1.5, fontface=3) +
  
  scale_color_manual(values=brewer.pal(11, "PiYG")[c(11, 3, 2, 1)]) +
  scale_y_continuous(limits=c(NA, 15)) +
  stat_pvalue_manual(data = pairwise.comparisons, y.position = 14, step.increase = 0.1, tip.length = 0, label = "p.signif", hide.ns = TRUE) +
  
  labs(x="ILF2 promoter mutation and Gain(1q) status", y="ILF2 Gene expression\n(VST normalized value)") +
  
  theme_bw(base_size = 6, base_line_size = 0.25, base_rect_size = 0.25) +
  theme(panel.grid = element_blank(), legend.position = "none")
ggsave("../figures/commpass_ilf2_any_mutant_gain1q.pdf", plot=boxplot.t.test.pairwise.ILF2.Gain1q, width = 2, height = 1.5) 
 