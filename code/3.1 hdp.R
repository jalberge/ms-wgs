setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
# 

sigs.long <- sigs |>
  tibble::rownames_to_column("SBS") |>
  pivot_longer(-SBS, names_to = "Tumor_Sample_Barcode", values_to = "SBS_Proportion")

# maf sbs:
maf.participants <- sort(unique(maf.sbs$Participant_ID))
maf.samples <- sort(unique(maf.sbs$SampleID))

mutation.burden <- maf.sbs %>% group_by(Tumor_Sample_Barcode) %>% count()
mutation.burden %>% write_tsv("../data/mut_burden_sample.txt")

sigs.df <- sigs.long

sigs.df <- sigs.df %>% inner_join(mutation.burden, by=c("Tumor_Sample_Barcode"))

# TODO LOOK AT THIS
# sensitivity=mean(unlist(lapply(rpois(n=100000,lambda=mean_cov),function(x) pbinom(q=2,size=x,p=VAF,lower.tail = F))))


# Bring participant level metadata to signatures sbs maf ------------------

ref.table <- clinical.and.terra.ref.pairs.samples %>% 
  select(Participant_ID, Reference_Pair, CTF_ID, HDP_Participant_ID, HDP_Reference_Pair)

maf.participants <- sort(unique(maf.sbs$Participant_ID))
matrix.participants <- sort(unique(meta$`Sample ID`))
lt <- list_to_matrix(list(MAF=maf.participants, META=matrix.participants))
m1 <-  make_comb_mat(lt)
UpSet(m1)

missing.in.maf <- as.data.frame(lt) |> filter(META==1 & MAF==0) |> rownames()
missing.in.matrix <- as.data.frame(lt) |> filter(META==0 & MAF==1) |> rownames()

# Rename meta participants with CTF ID when is CTF in SBS to make things easier

# If participant matches between MAF and META
match.meta <- meta %>% filter(`Sample ID` %in% maf.participants) # 940
nomatch.meta <- meta %>% filter(`Sample ID` %nin% maf.participants) # 92
nomatch.maf <- maf.participants[maf.participants %nin% meta$`Sample ID`] # 18
# If would match with CTF
ctf.matching <- ref.table %>% filter(!is.na(CTF_ID)) %>% filter(CTF_ID%in%nomatch.maf)
ctf.match.rescue <- data.frame(CTF_ID=nomatch.maf) %>% left_join(ctf.matching, by=c("CTF_ID")) %>% filter(!is.na(Participant_ID))
ctf.match.norescue <- data.frame(CTF_ID=nomatch.maf) %>% left_join(ctf.matching, by=c("CTF_ID")) %>% filter(is.na(Participant_ID))
rematch.ctf.meta <- nomatch.meta %>% 
  left_join(ctf.match.rescue, by=c("Sample ID"="Participant_ID")) %>% 
  filter(!is.na(CTF_ID)) %>% select(-`Sample ID`)%>%
  rename(`Sample ID`=CTF_ID) # LOGIC IS HERE

match.match.ctf.meta <- rbind(match.meta, rematch.ctf.meta %>% select(-Reference_Pair))

all(match.match.ctf.meta$`Sample ID` %in% maf.participants)
#TRUE
clinical.and.terra.ref.pairs.samples$Participant_ID[clinical.and.terra.ref.pairs.samples$Participant_ID %nin% match.match.ctf.meta$`Sample ID`]
# the only problem is that 13 patients are identified as CTF and SU2C.
# Honestly we should switch to PANGEA IDs now
# match match ctf meta contains participant IDs with mutations in MAF
# plus the problem comes from the meta file for matrix. which is not great.

# entered as Participant_ID
match.match.ctf.meta %>% write_tsv("../data/_MM_Sigs_HDP/20231206_HRD_ANNOT_MMRF_SU2C_WTC_CTF.tsv")

table(missing.in.maf %in%ref.table$Participant_ID)
table(missing.in.maf %in%ref.table$Reference_Pair)
table(missing.in.matrix %in%ref.table$CTF_ID)

missing.in.maf[missing.in.maf%in%ref.table$Participant_ID]
missing.in.maf[missing.in.maf%in%ref.table$Reference_Pair]
# missing CSL193, pM10993, pM11003, pM4175 => Are CLL etc.

# merge with --------------------------------------------------------------

maf.sbs <- maf.sbs |> left_join(match.match.ctf.meta, by=c("Participant_ID"="Sample ID"))

maf.sbs <- maf.sbs |> 
  mutate(Chromosome=factor(Chromosome, levels=c(1:22, "X", "Y"))) |>
  mutate(clonal=factor(clonal, levels=c(0, 1)))

# represent the data
quick.absolute.ccf <- function(maf, sample="CSL141") {
  maf |> 
    filter(Participant_ID==sample) |>
    ggplot(aes(Start_position)) +
    geom_pointrange(aes(y=multiplicity*ccf_hat, 
                        ymin=multiplicity*ccf_CI95_low, 
                        ymax=multiplicity*ccf_CI95_high,
                        color=clonal), alpha=0.6, size=0.25) +
    scale_color_manual(values=c("0"="lightgrey", "1"="#458B74")) +
    facet_grid(Tumor_Sample_Barcode~Chromosome, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position = "top") +
    labs(x="", title = paste0(sample, " ", "Multiplicity / clonality check"))
}


# Export pdf --------------------------------------------------------------

samples <- sort(match.match.ctf.meta$`Sample ID`)[1:5]

pdf("../figures/multiplicity_check.pdf")
for(sample in samples)
plot(quick.absolute.ccf(maf.sbs, sample))
dev.off()

# Check ratio multiplicity 1 vs 2 clonal mut ------------------------------


## Function def ------------------------------------------------------------

extract.hrd.counts <- function(maf, sample=NA, chromosomes=c(3, 5, 7, 9, 11, 15, 19, 21), min.mutations.per.multiplicity=10) {
  if (!is.na(sample)) maf <- maf |> 
      filter(Participant_ID==sample)
  maf |>
    filter(Chromosome %in% chromosomes) |>
    filter(clonal==1) |>
    group_by(Chromosome, multiplicity) |>
    count() |>
    group_by(Chromosome) |>
    filter(all(c(1, 2) %in% unique(multiplicity))) |>
    filter(n[multiplicity==2] > min.mutations.per.multiplicity & n[multiplicity==1] > min.mutations.per.multiplicity) |>
    filter(multiplicity %in% c(1,2))
}

compute.overall.ratio <-  function(counts) {
  counts |>
  ungroup() |>
    summarise(x=sum(n[multiplicity==2]), binom_test(x=sum(n[multiplicity==2]), n=sum(n[multiplicity<=2]))) |>
    mutate(Chromosome="HRD")
}

compute.chromosome.ratios <- function(counts) {
  counts |> 
    group_by(Chromosome) |> 
    summarise(x=n[multiplicity==2], binom_test(x=n[multiplicity==2], n=sum(n[multiplicity<=2])))
}

compute.ratios <- function(counts) {
  if(is.na(counts) | nrow(counts)<1) break
  r1 <- compute.overall.ratio(counts)
  r2 <- compute.chromosome.ratios(counts)
  bind_rows(r1, r2)
}

# plot.ratios <- function(ratios, Chromosomes=c(3, 5, 7, 9, 11, 15, 19, 21, "HRD")){

plot.ratios.hrd <- function(ratios){
  # select features
  ratios <- ratios |> filter(Chromosome == "HRD")
  # stage level statistics
  stats <- ratios |> 
    group_by(Stage) |> 
    summarise(mean=weighted.mean(x=estimate, w=1/(conf.high-conf.low)))
  #plot
  ggplot(ratios, aes(y=estimate, 
                     ymin=conf.low, 
                     ymax=conf.high, 
                     x=Participant_ID, 
                     color=Stage)) +
    geom_pointrange(size=0.25, shape=21) +
    geom_hline(data = stats, aes(yintercept = mean), color="black") +
    # scale_x_continuous(limits=c(0, 1)) +
    theme_pubclean() +
    facet_grid(~Stage, scales = "free_x") +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    labs(y="Fraction Clonal Mutations Multiplicity 2", title = "Clock-like estimate of acquisition of hyperdiploidy")
}

refactor.patients.hrd <- function(x) {
  x |>
    filter(Chromosome == "HRD") |>
    mutate(Stage=ifelse(Stage=="NDMM", "MM", Stage)) |>
    mutate(Participant_ID=fct_reorder(Participant_ID, estimate)) |>
    mutate(Stage=factor(Stage, levels=c("MGUS", "LRSMM", "IRSMM", "HRSMM", "MM")))
}


## Test run -----------------------------------------------

# missing ctf021, ctf034 ctf060 pM11077

hrd.patients <- match.match.ctf.meta |>
  filter(Cohort%in%c("SU2C", "CTC", "OBEN", "WTC", "MMRF")) |>
  filter(HRDTx=="HRD") |>
  # filter(`Sample ID` %in% maf.sbs$Tumor_Sample_Barcode) |> # FIXME here
  pull(`Sample ID`)

test.maf <- maf.sbs |> 
  filter(Participant_ID %in% hrd.patients)
  # filter(Tumor_Sample_Barcode==Tumor_Sample_Barcode_combined) # FOR NOW KEEP ONLY SINGLE SAMPLE PATIENTS

# hrd.patients <- hrd.patients[hrd.patients %in% test.maf$Tumor_Sample_Barcode]

hrd.counts <- test.maf |>
  split(test.maf$Participant_ID) |>
  map(\(df) extract.hrd.counts(df, min.mutations.per.multiplicity = 10)) |>
  list_rbind(names_to = "Participant_ID") |> 
  left_join(match.match.ctf.meta, by=c("Participant_ID"="Sample ID"))

hrd.ratios <- hrd.counts |>
  split(hrd.counts$Participant_ID) |>
  map(\(cts) compute.ratios(cts)) |>
  list_rbind(names_to = "Participant_ID") |>
  left_join(match.match.ctf.meta, by=c("Participant_ID"="Sample ID"))

ratios.hrd <- plot.ratios.hrd(refactor.patients.hrd(hrd.ratios |> filter(Cohort!="MMRF")))
ratios.hrd

ggsave("../figures/multiplicity_estimate_hyperdiploid.png", ratios.hrd, width = 4, height = 2, scale = 2)

# match annot for sig weight ----------------------------------------------

sigs.df <- sigs.df %>%
  mutate(Sample=case_when(startsWith(Tumor_Sample_Barcode, "MMRF_") ~ str_extract(Tumor_Sample_Barcode, "MMRF_[0-9]{4}"), TRUE ~ Tumor_Sample_Barcode))

missing.samples.in.meta.matrix <- unique(sigs.df$Sample[ sigs.df$Sample %nin% meta.sigs$`Sample ID` ])
# convert CTF to CAD ID

# what is the point of this table
ref.table <- clinical.and.terra.ref.pairs.samples %>% 
  select(Participant_ID, Reference_Pair)
  # filter(Reference_Pair %nin% meta.sigs$`Sample ID`)

sigs.df <- sigs.df %>%
  left_join(ref.table, by=c("Sample"="Reference_Pair"))
  # mutate(Sample = case_when(!is.na(Participant_ID)~Participant_ID, TRUE~Sample))
  # select(-Participant_ID)
# done!
sigs.df$Sample[ sigs.df$Sample %nin% meta.sigs$`Sample ID` ]

sigs.df <- sigs.df %>% filter(Sample %nin% c("PRM_NEG4622", "pM10993")) # hihi mgip signature

meta.sigs$`Sample ID`[ meta.sigs$`Sample ID` %nin% sigs.df$Sample ]

# not sure what happened here ...
#> meta.sigs$`Sample ID`[ meta.sigs$`Sample ID` %nin% sigs.df$Sample ]
# [1] "MMRF_1286" "MMRF_2330" "MMRF_1342" "MMRF_1510" "MMRF_2272" "MMRF_1957" "MMRF_2433" "MMRF_1641" "MMRF_2595"
# [10] "MMRF_1646" "MMRF_1715" "MMRF_1726" "MMRF_1847" "MMRF_2507" "MMRF_2089" "MMRF_2659" "MMRF_2722" "MMRF_2833"

metasigs <- sigs.df %>%
  inner_join(meta.sigs, by=c("Sample"="Sample ID")) %>%
  mutate(across(where(is.double), as.numeric))

refit.metasigs <- metasigs %>% 
  mutate(sum_SBS_mean=rowSums(select(., matches("^SBS.*_mean$"))),
         sum_SBS_lower95=rowSums(select(., matches("^SBS.*_lower95$"))),
         sum_SBS_upper95=rowSums(select(., matches("^SBS.*_upper95$")))) %>%
  mutate(adjusted_n=round(sum_SBS_mean*n)) %>%
  mutate(across(matches("^SBS.*_mean$"), ~.x/sum_SBS_mean)) %>%
  mutate(across(matches("^SBS.*_lower95$"), ~.x/sum_SBS_lower95)) %>%
  mutate(across(matches("^SBS.*_upper95$"), ~.x/sum_SBS_upper95))

flat.refit <- refit.metasigs %>% 
  pivot_longer(cols=matches("^(SBS|N).*(_mean|_lower95|_upper95)"), 
               names_to = c("SBS", "measure"), 
               names_pattern = "(.*)_(.*)") %>% 
  pivot_wider(names_from = measure, values_from = value)

flat.refit.sum <- flat.refit %>%
  mutate(across(c(mean, lower95, upper95), ~.x*adjusted_n))

flat.refit.sum %>% write_tsv("../data/_MM_Sigs_HDP/Flat_Refit_Sum.tsv")





# sigs --------------------------------------------------------------------

flat.refit.sum <- read_tsv("../data/_MM_Sigs_HDP/Flat_Refit_Sum.tsv")

# quick Total N mutations compute

flat.refit.sum %>% filter(!startsWith(Sample, "MMRF")) %>% summary(adjusted_n)

colnames(flat.refit)

flat.refit.sum %>% 
  filter(startsWith(SBS, "SBS")) %>%
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(Age=as.numeric(Age)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  ggplot(aes(Age, mean, color=Stage)) + 
  geom_point(shape=21, alpha=0.5) + 
  geom_smooth(method = "lm") +
  facet_wrap(~SBS, scales = "free")

flat.refit.sum %>% 
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(Age=as.numeric(Age)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  group_by(Sample, Stage, Age) %>%
  summarise(`SBS_5_40`=mean[SBS=="SBS5"]+mean[SBS=="SBS40"]) %>% 
  filter(!is.na(SBS_5_40)) %>%
  ggplot(aes(Age, SBS_5_40, color=Stage)) + 
  geom_point(shape=21, alpha=0.5) + 
  geom_smooth(method = "lm") +
  facet_wrap(~Stage, scales = "free")

pt.order.apobec <- flat.refit.sum %>% filter(SBS%in%c("SBS2", "SBS13")) %>% group_by(Sample) %>% summarise(sum=sum(mean)) %>% arrange(sum) %>% pull(Sample)

flat.refit.sum %>% 
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(Age=as.numeric(Age)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  group_by(Sample, Stage, Age, HRDTx) %>%
  summarise(`SBS_APOBEC`=mean[SBS=="SBS2"]+mean[SBS=="SBS13"]) %>% 
  filter(!is.na(SBS_APOBEC)) %>%
  ggplot(aes(factor(Sample, pt.order.apobec), SBS_APOBEC, fill=Stage)) + 
  geom_col() + 
  facet_wrap(~HRDTx, scales = "free")

flat.refit %>% 
  filter(startsWith(SBS, "SBS")) %>%
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(Age=as.numeric(Age)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  ggplot(aes(Age, mean, color=Stage)) + 
  geom_point(shape=21, alpha=0.5) + 
  geom_smooth(method = "lm") +
  facet_wrap(~SBS, scales = "free")

# Plot --------------------------------------------------------------------

palette <- c("SBS1"="#A6CEE3", "SBS2"="#FB9A99", "SBS5"="#1F78B4", "SBS8"="#B2DF8A", "SBS9"="#33A02C",
             "SBS13"="#E31A1C", "SBS16"="#FDBF6F", "SBS17a"="#CAB2D6", "SBS17b"="#6A3D9A", "SBS18"="#FF7F00",
             "SBS19"="#FFFF99", "SBS40"="#6BAED6")
# blues
# c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", 
#   "#2171B5", "#084594")

# pt.order <- flat.refit %>% filter(SBS=="SBS9") %>% arrange(mean) %>% pull(Sample)
pt.order <- flat.refit.sum %>% group_by(Sample) %>% slice_max(adjusted_n, with_ties = FALSE) %>% arrange(-adjusted_n) %>% pull(Sample)

flat.refit %>% 
  filter(startsWith(SBS, "SBS")) %>%
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  ggplot(aes(factor(Sample, levels=pt.order), mean, fill=fct_reorder(SBS, mean))) + 
  geom_col(position="stack") +
  scale_fill_manual(values=palette) +
  scale_y_continuous(expand = expansion(0, 0)) +
  facet_grid(~Stage, scales = "free_x") +
  labs(x="Participant", y="SBS mean weight", fill="SBS") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(),
        panel.grid = element_blank())

flat.refit.sum %>% 
  filter(startsWith(SBS, "SBS")) %>%
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  mutate(SBS=factor(SBS, levels=paste0("SBS", c(1, 2, 5, 8, 9, 13, 16, "17a", "17b", 18, 19, 40)))) %>%
  mutate(Stage=factor(Stage, levels=c("MGUS", "SMM", "MM"))) %>%
  ggplot(aes(factor(Sample, levels=pt.order), mean, fill=fct_reorder(SBS, mean))) + 
  geom_col(position="stack") +
  scale_fill_manual(values=palette) +
  scale_y_continuous(expand = expansion(0, 0)) +
  facet_grid(~Stage, scales = "free_x") +
  labs(x="Participant", y="SBS mean weight", fill="SBS") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(),
        panel.grid = element_blank())


flat.refit %>% 
  filter(SBS=="SBS9") %>%
  mutate(Stage=case_when(Stage%in%c("LRSMM", "IRSMM", "HRSMM")~"SMM", 
                         Stage=="NDMM"~"MM",
                         TRUE~Stage)) %>%
  group_by(Stage) %>%
  count

ggplot(metasigs, aes(Stage, SBS1_mean)) + geom_boxplot()
ggplot(refit.metasigs, aes(Stage, SBS1_mean+SBS5_mean)) + geom_boxplot()

flat.refit %>%
  filter(startsWith(SBS, "SBS") & Cohort != "MMRF") %>%
  ggplot(aes(factor(Stage, levels=c("MGUS", "LRSMM", "IRSMM", "SMM", "HRSMM", "NDMM", "MM")), mean, color=Stage)) + 
  geom_boxplot() + facet_wrap(~SBS, scales = "free")

flat.refit %>%
  filter(startsWith(SBS, "SBS")) %>%
  ggplot(aes(HRDTx, mean, color=HRDTx)) + geom_boxplot() + geom_jitter(alpha=0.2) + facet_wrap(~SBS)

flat.refit %>%
  filter(startsWith(SBS, "SBS")) %>%
  ggplot(aes(Gender, mean, color=Gender)) + geom_boxplot() + geom_jitter(alpha=0.2) + facet_wrap(~SBS)

flat.refit %>%
  filter(startsWith(SBS, "SBS")) %>%
  ggplot(aes(Cohort, mean, color=Cohort)) + geom_boxplot() + geom_jitter(alpha=0.2) + facet_wrap(~SBS) +
  stat_compare_means()

ggplot(metasigs, aes(Stage, N13_mean)) + geom_boxplot()


# signature in ig regions -------------------------------------------------

most.likely.mutation <- maf.sbs %>%
  filter(Chromosome==14 & Start_position >= 106320000 & Start_position <= 106335000) %>%
  mutate(sum_SBS_mean=rowSums(select(., matches("^SBS[0-9ab]{1,}")))) %>%
  mutate(across(matches("^SBS[0-9ab]{1,}"), ~.x/sum_SBS_mean)) %>%
  pivot_longer(cols=matches("^(SBS|N)[0-9ab]{1,}"), 
               names_to = c("SBS"), 
               values_to = "weight") %>% 
  filter(startsWith(SBS, "SBS")) %>% 
  group_by(Tumor_Sample_Barcode, Chromosome, Start_position) %>%
  slice_max(n=1, weight)

IGH.fractions <- most.likely.mutation %>%
  group_by(Tumor_Sample_Barcode, SBS) %>%
  count() %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(Frac=n/sum(n))

meta.igh <- IGH.fractions %>% inner_join(flat.refit, by=c("Tumor_Sample_Barcode"="Sample", "SBS"))

ggplot(meta.igh, aes(HRDTx, Frac, color=HRDTx)) + geom_boxplot() + geom_jitter(alpha=0.2) + facet_wrap(~SBS)

