setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# code to list, classify, and compare non-coding candidate drivers from DIG analysis
# this includes adding normal B cells, spliting categories of translocations, etc.

library(data.table)
library(tidyverse)

source("0_annotate_samples.R")


# Load DIG ref ------------------------------------------------------------

pcawg.bed.hg19 <- read_table("../data/_DIG/MM/annotations/grch37.PCAWG_noncoding.bed", 
                             col_names = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts" ))
nc.pcawg.bed.hg19 <- pcawg.bed.hg19 |> filter(!str_detect(name, "gc19_pc.cds"))


# Load DIG results --------------------------------------------------------

nc.dig.results <- read_table("../data/_DIG/MM/unfiltered_dig_results/su2c-commpass/su2c-179-commpass-856-hg19-dig.pcawg_noncoding.results.txt")
nc.dig.results <- nc.dig.results |> 
  filter(!str_detect(ELT, "gc19_pc.cds"))  |>
  mutate(q=p.adjust(PVAL_MUT_BURDEN, method="bonferroni")) |>
  arrange(q) |>
  filter(q<0.05)

all(nc.dig.results$ELT %in% nc.pcawg.bed.hg19$name)


# Add coordinates to DIG results ------------------------------------------

nc.dig.results.bed <- nc.dig.results |> left_join(nc.pcawg.bed.hg19, by=c("ELT"="name"))

nc.dig.results.bed <- nc.dig.results.bed |> separate(ELT, sep="::", into=c("elt_group", "elt_source", "elt_gene_name", "elt_identifier"), remove = FALSE)

nc.dig.results.bed |> write_tsv("../data/nc.dig.results.bed.tsv")

# Load MAF ----------------------------------------------------------------

debug=FALSE
if(debug) { n_max=100000 } else { n_max =Inf }

maf.path <- "../data/_MM_Sigs_HDP/20240304_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level_combined_signatures_annotated.tsv.gz"
maf.sbs <- read_tsv(maf.path, num_threads = 7, n_max = n_max)


# Extract maf for one element ---------------------------------------------

test.nc.dig.results <- nc.dig.results.bed[1,]

filter.maf <- function(bed, maf) {
  # print(bed)
  # message(bed["chr"], ":", bed["start"], "-", bed["end"])
  # message(class(bed))
  print(pull(bed, chr))
  maf <- maf |>
    filter(Chromosome==pull(bed, chr) & Start_position>=pull(bed, start) & End_position<=pull(bed, end))
  maf |>
    mutate(elt_chr=pull(bed, chr), elt_start=pull(bed, start), elt_end=pull(bed, end)) |>
    left_join(bed, by=c("elt_chr"="chr", "elt_start"="start", "elt_end"="end"))
}

library(tictoc)

#0.19 for 1 element
# 2 sec for all elements on 100;000 mutations
# 8 sec for all elements on 2E6 mutations

tic()
# nc.maf <- apply(nc.dig.results.bed[1,], 1, filter.maf, maf=maf.sbs)
nc.maf <- pmap(nc.dig.results.bed, function(...){x=tibble(...);return(filter.maf(x, maf.sbs))}) |> list_rbind()
toc()


# Add annotations for patients --------------------------------------------

meta <- read_tsv("../data/FINAL_meta_columns_matrix.tsv")

meta.hdp <- meta %>% 
  filter(Cohort %nin% c("Bustoros")) %>%
  left_join(clinical.and.terra.ref.pairs.samples %>% select(HDP_Participant_ID, HDP_Reference_Pair, Participant_ID), by=c("Sample ID" = "Participant_ID")) %>%
  mutate(HDP_Participant_ID=ifelse(is.na(HDP_Participant_ID), `Sample ID`, HDP_Participant_ID))

nc.maf <- nc.maf |> inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID"))
maf.sbs <- maf.sbs |> inner_join(meta.hdp, by=c("Participant_ID"="HDP_Participant_ID"))

nc.maf |> filter(elt_gene_name=="WDR74") |> write_tsv("../data/review_wdr74_nc.tsv")

# Signature summary per element -------------------------------------------

sbs.palette = 
  c("N9"="#333333", "N12"="#5c5c5c", "N16"="#858585", "N18"="#adadad", "N19"="#d6d6d6", "N20"="#ffffff",
    "N13"="#C7E9C0", "N15"="#A1D99B", "SBS8"="#74C476", "SBS16"="#31A354", "SBS18"="#006D2C",
    "SBS1"="#FA8072", "SBS5"="#B22222", "SBS40"="#B22222", "SBS2"="#FFA500", "SBS13"="#EEEE00", 
    "SBS17a"="#EEAEEE", "SBS17b"="#B452CD", "SBS9"="#BFEFFF", "N17"="#00BFFF", "SBS85"="#1874CD")
sbs.order <- c("N9", "N12", "N16", "N18", "N19", "N20", "N13", "N15", "SBS8", 
               "SBS16", "SBS18", "SBS1", "SBS5", "SBS40", "SBS2", "SBS13", "SBS17a", 
               "SBS17b", "SBS9", "N17", "SBS85")
# signature
nc.maf.sbs <- nc.maf |> 
  group_by(ELT, elt_chr, elt_start, elt_end, elt_group, elt_source, elt_gene_name, elt_identifier) |>
  summarise_at(vars(matches("^(N|SBS)[0-9]{1,}|^SBS[0-9]{1,}")), mean) |>
  pivot_longer(-c(ELT, elt_chr, elt_start, elt_end, elt_group, elt_source, elt_gene_name, elt_identifier), names_to = "SBS", values_to = "weight") |>
  mutate(SBS=factor(SBS, levels=sbs.order))

nc.maf.sbs |> write_tsv("../data/dig_input_sbs_20230314.tsv")
nc.maf.sbs <-  read_tsv("../data/dig_input_sbs_20230314.tsv")

axis.labels <- nc.maf.sbs |> 
  select(ELT, elt_chr, elt_start, elt_end, elt_group, elt_source, elt_gene_name, elt_identifier) |> 
  distinct() |> 
  mutate(label=paste0(elt_gene_name, "(", elt_identifier, ") ", elt_chr, ":", elt_start, "-", elt_end))
# paste0(elt_gene_name, "(", elt_identifier, ")")

sbs.nc <- ggplot(nc.maf.sbs, aes(weight, ELT, fill=SBS)) +
  
  geom_bar(stat="identity", position = "fill") +
  
  scale_fill_manual(values=sbs.palette) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_discrete(breaks=axis.labels$ELT, labels=axis.labels$label) +
  
  theme_minimal_vgrid(font_size = 6) +
  theme(panel.grid = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

# subset elts of interest in maf

ggsave("../figures/sbs.nc.dig.png", width = 3, height = 4, scale=2)
ggsave("../figures/sbs.nc.dig.pdf", width = 3, height = 4, scale=2)
# make sure sigs are there
# make sure ccf is there, too



# Split RNA / regulatory --------------------------------------------------

# small rna and lnc <- 

rna.groups <- c("smallrna.ncrna", "mirna_pre", "mirna_mat", "mirna.prom", "lncrna.promCore", "lncrna.ncrna")
regulatory <- c("gc19_pc.promCore", "gc19_pc.5utr", "gc19_pc.3utr")
# drop enhancers as only Ig + BMP6 already in promoter

length(unique(nc.dig.results.bed$ELT))
length(unique(nc.maf.sbs$ELT))

unique(nc.dig.results.bed$ELT)[unique(nc.dig.results.bed$ELT) %nin%unique(nc.maf.sbs$ELT) ]

# pick promoters and UTRs and enhancers
nc.maf.sbs <- nc.maf.sbs |>
  filter(elt_group %in% regulatory) |>  
  filter(elt_gene_name!="IGLL5") |>
  filter(!(elt_chr==14 & elt_start > 105E6))


# # for artefacts
# ref.depth <- maf.sbs |>
#   filter(ccf_hat>0.9) |>
#   head(n = 100000) |>
#   mutate(depth=t_alt_count+t_ref_count)
# 
# ggplot(ref.depth, aes(depth)) + geom_histogram()
# 
# nc.maf.ref <- nc.maf |> 
#   filter(ELT=="gc19_pc.promCore::gencode::IGLL5::ENSG00000254709.2") |>
#   mutate(depth=t_alt_count+t_ref_count)
# 
# nc.maf |> 
#   mutate(depth=t_alt_count+t_ref_count) |>
#   group_by(ELT) |>
#   summarize(p=wilcox.test(depth, nc.maf.ref$depth)$p.value) |>
#   View()



# Show depth of sequencing ------------------------------------------------

nc.maf$depth <- nc.maf$t_alt_count+nc.maf$t_ref_count

nc.maf.depth <- nc.maf |> 
  filter(!startsWith(Tumor_Sample_Barcode, "MMRF")) |>
  group_by(ELT, Participant_ID) |>
  slice_sample(n=1) |>
  group_by(ELT) |> 
  mutate(depth=t_alt_count+t_ref_count) |> 
  # mutate(ELT=fct_reorder(ELT, depth, mean, .na_rm=TRUE)) |>
  ggplot(aes(depth, fct_reorder(ELT, depth, mean, .na_rm=TRUE))) + 
  geom_boxplot() +
  geom_vline(xintercept = 60) +
  scale_x_continuous()
nc.maf.depth

ggsave("../figures/nc.maf.depth.pdf", width = 6, height = 4, scale=2)

# # test depth of sequencing
# nc.maf |>
#   group_by(ELT) |>
#   summarize(p=ks.test(depth, "ppois", 60, alternative = "l")$p.value) |>
#   View()

# list of artefacts
# after reviewing mutations in IGV
artefact.elements <- c("gc19_pc.3utr::gencode::PPP1CA::ENSG00000172531.10",
  "gc19_pc.5utr::gencode::MYO1E::ENSG00000157483.4", 
  "gc19_pc.5utr::gencode::NFS1::ENSG00000244005.8",
  "gc19_pc.5utr::gencode::WDR74::ENSG00000133316.11", 
  "gc19_pc.promCore::gencode::FAM81A::ENSG00000157470.7", 
  "gc19_pc.promCore::gencode::GRIK5::ENSG00000105737.5",
  "gc19_pc.promCore::gencode::WDR74::ENSG00000133316.11")

# filter out artefacts
nc.maf.sbs <- nc.maf.sbs |> 
  filter(ELT %nin% artefact.elements)

# SBS non coding ----------------------------------------------------------


# reorder by APOBEC
order.apobec <- nc.maf.sbs |> 
  # filter(SBS=="SBS2") |>
  filter(SBS %in% c("SBS9", "N17", "SBS85")) |>
  group_by(ELT) |>
  summarise(weight=sum(weight)) |>
  arrange(-weight) |>
  pull(ELT)

# add coordinates to labels for now
axis.labels <- nc.maf.sbs |> 
  select(ELT, elt_chr, elt_start, elt_end, elt_group, elt_source, elt_gene_name, elt_identifier) |> 
  distinct() |> 
  mutate(label=paste0(elt_gene_name, "(", elt_identifier, ") ", elt_chr, ":", elt_start, "-", elt_end))

# axis.labels$label <- factor(axis.labels$label, levels=axis.labels$label)

axis.labels <- axis.labels |>
  mutate(final.label= paste0( elt_gene_name, " ", str_remove(elt_group, "gc19_pc."))) |>
  mutate(final.label=str_replace(final.label, "utr", "'UTR"),) |>
  mutate(final.label=str_replace(final.label, "promCore", "Prom."),) |>
  select(ELT, final.label)

# reorder rows
nc.maf.final <- nc.maf |> filter(ELT %in%nc.maf.sbs$ELT) |> left_join(axis.labels)
nc.maf.final |> write_tsv("../data/export_maf_sbs_noncoding_regulatory.maf")
nc.maf.final <- read_tsv("../data/export_maf_sbs_noncoding_regulatory.maf")

## fig -------------

plot.sbs.nc.reg <- ggplot(nc.maf.sbs, aes(weight, factor(ELT, levels=order.apobec), fill=SBS)) +
  
  geom_bar(stat="identity", position = "fill") +
  
  scale_fill_manual(values=sbs.palette) +
  
  scale_x_continuous(limits=c(0, 1), labels=scales::percent) +
  scale_y_discrete(breaks=axis.labels$ELT, labels=paste0(axis.labels$final.label)) +
  
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), 
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank()
        )
  
  # facet_grid(elt_group ~ ., scales="free_y", space = "free_y")

plot.sbs.nc.reg

ggsave("../figures/sbs.nc.dig.reg.only.png", width = 3, height = 4, scale=2)
ggsave("../figures/sbs.nc.dig.reg.only.pdf", width = 3, height = 4, scale=2)


# Number of participants non coding ---------------------------------------

N_participants = maf.sbs |> filter(!startsWith(Participant_ID, "MMRF")) |> select(Participant_ID) |> distinct() |> ungroup() |> count()

N_per_class <- maf.sbs |> filter(!startsWith(Participant_ID, "MMRF")) |> select(Participant_ID, IMWG) |> distinct() |> group_by(IMWG) |> count()

nc.frequencies <- nc.maf.final |> 
  filter(!startsWith(Participant_ID, "MMRF")) |>
  group_by(ELT) |>
  summarise(n=pull(N_participants, n),
            N_mut_participants=length(unique(Participant_ID)),
            Frac_mut_participants=N_mut_participants/pull(N_participants, n))
nc.frequencies$IMWG <- "Overall"

nc.frequencies.per.group <- nc.maf.final |> 
  filter(!startsWith(Participant_ID, "MMRF")) |>
  group_by(ELT, IMWG) |>
  summarise(N_mut_participants=length(unique(Participant_ID))) |>
  full_join(N_per_class) |>
  mutate(Frac_mut_participants=N_mut_participants/n)

nc.frequencies.global.and.per.group <- bind_rows(nc.frequencies, nc.frequencies.per.group)

nc.frequencies.global.and.per.group <- nc.frequencies.global.and.per.group |>
  mutate(IMWG = factor(IMWG, levels = c("Overall", "Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"))) |>
  mutate(Freq=cut(Frac_mut_participants, breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)))


## fig Overall -------------

plot.nc.overall.frequencies <- ggplot(nc.frequencies.global.and.per.group |> filter(IMWG=="Overall"), 
                                      aes(IMWG, factor(ELT, levels=order.apobec), 
                                          fill=Frac_mut_participants)) +
  geom_tile() +
  
  scale_y_discrete(breaks=axis.labels$ELT, labels=paste0(axis.labels$final.label)) +
  
  scale_fill_distiller(direction = 1, palette="YlGnBu", na.value = "white", limits=c(0, 1), labels=scales::percent) +
  
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_bw(base_size = 6) +
  
  labs(x="", y="") +
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.title.y = element_blank()
        ) +
  coord_equal()

## fig Per subgroup ------------

plot.nc.group.frequencies <- 
ggplot(nc.frequencies.global.and.per.group |> filter(IMWG!="Overall"), aes(factor(IMWG, 
                                                                                  levels = c("Unclassified", "Hyperdiploid", "Cyclin D", "MMSET", "MAF"), 
                                                                                  labels = c("Unc.", "HRD", "CCND", "MMSET", "MAF")), 
                                                                           factor(ELT, levels=order.apobec), fill=Frac_mut_participants)) +
  
  # scale_fill_gradient2(low = "white", mid="royalblue", high = "darkblue", na.value = "white", midpoint = 0.5) +
  
  geom_tile() +
  
  scale_y_discrete(breaks=axis.labels$ELT, labels=paste0(axis.labels$final.label)) +
  
  # scale_fill_fermenter(direction = 1, palette="YlOrRd") +
  scale_fill_distiller(direction = 1, palette="YlGnBu", na.value = "white", limits=c(0, 1), labels=scales::percent) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        
        legend.position = "top",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
        ) +
  coord_equal()


# Add qvalues -------------------------------------------------------------

nc.qvalues <- nc.dig.results.bed |>
  filter(ELT %in% axis.labels$ELT) |>
  mutate(`-log10(q-value)`=-log10(q))

plot.nc.qvalues <- ggplot(nc.qvalues, aes(`-log10(q-value)`, factor(ELT, levels=order.apobec))) +
  geom_bar(stat="identity", color="black", fill="white", width = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype=2, color="salmon", linewidth=0.5) +
  scale_y_discrete(breaks=axis.labels$ELT, labels=paste0(axis.labels$final.label)) +
  scale_x_continuous(trans="pseudo_log") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), 
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
  ) 

plot.nc.qvalues


# Add normal B cells ------------------------------------------------------

simple.b.hmach.maf <-  read_tsv("../data/20230228_try_vcf2maf_conv_normal_B_colonies_Machado.maf.gz")

# be fair compared to SBS and include only SNPs
simple.b.hmach.maf <- simple.b.hmach.maf |> filter(Variant_Type=="SNP")

table(simple.b.hmach.maf$Tumor_Seq_Allele2)
table(simple.b.hmach.maf$Reference_Allele)
table(simple.b.hmach.maf$Tumor_Seq_Allele2, simple.b.hmach.maf$Reference_Allele)

## counts mutants in bed ------------

nc.b.maf <- pmap(nc.dig.results.bed |> filter(ELT %in% axis.labels$ELT), function(...){x=tibble(...);return(filter.maf(x, simple.b.hmach.maf))}) |> list_rbind()

N_B_total <- simple.b.hmach.maf |> select(sample) |> distinct() |> count() |> pull(n)
N_B_memory_total <- simple.b.hmach.maf |> filter(Cell.type2=="Memory B") |> select(sample) |> distinct() |> count() |> pull(n)

N_B_per_class <- simple.b.hmach.maf |> select(sample, Cell.type2) |> distinct() |> group_by(Cell.type2) |> count()

nc.b.frequencies <- nc.b.maf |> 
  ungroup() |>
  full_join(axis.labels, by=c("ELT")) |>
  group_by(ELT, Cell.type2) |>
  summarise(N_mut_participants=length(unique(sample))) |>
  full_join(N_B_per_class) |>
  mutate(Frac_mut_participants=N_mut_participants/n)

# add missing values lost in joint
nc.b.frequencies <- nc.b.frequencies |> 
  mutate(Cell.type2 = replace_na(Cell.type2, "Memory B")) |>
  mutate(n = replace_na(n, 0)) |>
  mutate(Frac_mut_participants = replace_na(Frac_mut_participants, 0)) |>
  ungroup() |>
  complete(ELT, Cell.type2, fill = list(Frac_mut_participants = 0, n = 0))

nc.b.frequencies |> write_tsv("../data/nc_mut_28_drivers_utr_prom_freq_b_cells.tsv")

nc.b.frequencies |> filter(Cell.type2=="Memory B") |>
  mutate(Recurrent = N_mut_participants  > 1) |>
  group_by(Recurrent) |>
  count()

nc.b.frequencies |>
  filter(Cell.type2=="Memory B" & N_mut_participants  > 1)

# 17 elements non-coding

# A tibble: 17 × 5
# ELT                                                     Cell.type2 N_mut_participants     n Frac_mut_participants
# <chr>                                                   <chr>                   <int> <int>                 <dbl>
#   1 gc19_pc.5utr::gencode::BACH2::ENSG00000112182.10        Memory B                   25    74                0.338 
# 2 gc19_pc.5utr::gencode::BCL6::ENSG00000113916.13         Memory B                   43    74                0.581 
# 3 gc19_pc.5utr::gencode::BCL7A::ENSG00000110987.4         Memory B                   12    74                0.162 
# 4 gc19_pc.5utr::gencode::DNMT1::ENSG00000130816.10        Memory B                   15    74                0.203 
# 5 gc19_pc.5utr::gencode::LPP::ENSG00000145012.8           Memory B                   28    74                0.378 
# 6 gc19_pc.5utr::gencode::LRMP::ENSG00000118308.10         Memory B                    6    74                0.0811
# 7 gc19_pc.5utr::gencode::RFTN1::ENSG00000131378.9         Memory B                    3    74                0.0405
# 8 gc19_pc.5utr::gencode::SERPINA9::ENSG00000170054.10     Memory B                   12    74                0.162 
# 9 gc19_pc.promCore::gencode::ACTB::ENSG00000075624.9      Memory B                    5    74                0.0676
# 10 gc19_pc.promCore::gencode::BACH2::ENSG00000112182.10    Memory B                    4    74                0.0541
# 11 gc19_pc.promCore::gencode::BCL6::ENSG00000113916.13     Memory B                   43    74                0.581 
# 12 gc19_pc.promCore::gencode::IKZF3::ENSG00000161405.12    Memory B                    5    74                0.0676
# 13 gc19_pc.promCore::gencode::LPP::ENSG00000145012.8       Memory B                   32    74                0.432 
# 14 gc19_pc.promCore::gencode::LRMP::ENSG00000118308.10     Memory B                    6    74                0.0811
# 15 gc19_pc.promCore::gencode::PAX5::ENSG00000196092.8      Memory B                   13    74                0.176 
# 16 gc19_pc.promCore::gencode::RFTN1::ENSG00000131378.9     Memory B                    5    74                0.0676
# 17 gc19_pc.promCore::gencode::SERPINA9::ENSG00000170054.10 Memory B                    6    74                0.0811


nc.b.frequencies |>
  filter(Cell.type2=="Memory B" & N_mut_participants  <= 1)
# # A tibble: 11 × 5
# ELT                                                   Cell.type2 N_mut_participants     n Frac_mut_participants
# <chr>                                                 <chr>                   <int> <int>                 <dbl>
#   1 gc19_pc.3utr::gencode::BMP6::ENSG00000153162.8        Memory B                    1     0                0     
# 2 gc19_pc.3utr::gencode::BTG2::ENSG00000159388.5        Memory B                    1     0                0     
# 3 gc19_pc.3utr::gencode::SYBU::ENSG00000147642.12       Memory B                    1     0                0     
# 4 gc19_pc.3utr::gencode::TXNDC5::ENSG00000239264.4      Memory B                    1     0                0     
# 5 gc19_pc.5utr::gencode::CCND1::ENSG00000110092.3       Memory B                    1     0                0     
# 6 gc19_pc.5utr::gencode::DTX1::ENSG00000135144.3        Memory B                    1    74                0.0135
# 7 gc19_pc.5utr::gencode::EGR1::ENSG00000120738.7        Memory B                    1    74                0.0135
# 8 gc19_pc.promCore::gencode::CCND1::ENSG00000110092.3   Memory B                    1     0                0     
# 9 gc19_pc.promCore::gencode::DTX1::ENSG00000135144.3    Memory B                    1     0                0     
# 10 gc19_pc.promCore::gencode::ILF2::ENSG00000143621.12   Memory B                    1    74                0.0135
# 11 gc19_pc.promCore::gencode::STRIP1::ENSG00000143093.10 Memory B                    1     0                0     


## fig -----------------


nc.b.frequencies.plot <- 
  ggplot(nc.b.frequencies, aes(factor(Cell.type2, levels=c("Naive B", "Memory B")), factor(ELT, levels=order.apobec), fill=Frac_mut_participants)) +
  
  # scale_fill_gradient2(low = "white", mid="royalblue", high = "darkblue", na.value = "white", midpoint = 0.5) +
  
  geom_tile() +
  
  scale_y_discrete(breaks=axis.labels$ELT, labels=paste0(axis.labels$final.label)) +
  
  scale_fill_gradient2(low = "white", high = "darkblue", na.value = "white", limits=c(0, 1), labels=scales::percent) +
  # scale_fill_distiller(direction = 1, palette="Purples", na.value = "white", limits=c(0, 1), labels=scales::percent) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  
  labs(x="", y="", fill="% B cell mut.") +
  
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        
        legend.position = "top",
        # axis.line.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.title.y = element_blank(),
  ) +
  coord_equal()

nc.b.frequencies.plot

ggsave("../figures/B_cell_mutation_hotspots.png")

# Assemble figure ---------------------------------------------------------

plot.sbs.for.final <- plot.sbs.nc.reg +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) 

plot.nc.overall.frequencies.for.final <- plot.nc.overall.frequencies +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) 

plot.nc.group.frequencies.for.final <- plot.nc.group.frequencies +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position = "top",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank()) 


library(patchwork)

final.non.coding.matrix <- nc.b.frequencies.plot + 
  plot.nc.overall.frequencies.for.final + 
  plot.nc.group.frequencies.for.final + 
  plot.sbs.for.final + 
  plot.nc.qvalues +
  plot_layout(nrow = 1, widths = c(1, 1, 1, 1, 0.4))

final.non.coding.matrix

ggsave2("../figures/noncodingmatrix_patchwork.pdf", final.non.coding.matrix, width = 4, height = 4, scale=1.5)
