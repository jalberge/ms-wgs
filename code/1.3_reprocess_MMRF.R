setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

source("0_annotate_samples.R")

# README ------------------------------------------------------------------

# tx come from CatchTheFISH + SV_WORKFLOW
# CNV come from wgs workflow on hg38
# genome-wide mutation come from wgs workflow on hg38
# exome mutation come from rsp's analysis of exomes
# remove samples that were removed from Keats analysis (portal?) or MMRF IA20

# Clinical ----------------------------------------------------------------
# here sequencing number starts at 1 (like MMRF_1234_1) but we rely on Baseline VJ annotation to extract correct sample date / time / screening
mmrf.pt <- read_tsv("../../../MMRF_IA20/CoMMpass_IA20_FlatFiles/MMRF_CoMMpass_IA20_PER_PATIENT.tsv")
mmrf.pt.visit <- read_tsv("../../../MMRF_IA20/CoMMpass_IA20_FlatFiles/MMRF_CoMMpass_IA20_PER_PATIENT_VISIT.tsv")
mmrf.pt.baseline.visit.with.seq <- mmrf.pt.visit %>% filter(!is.na(SPECTRUM_SEQ)) %>% filter(VJ_INTERVAL=="Baseline")
mmrf.pt.os <- read_tsv("../../../MMRF_IA20/CoMMpass_IA20_FlatFiles/MMRF_CoMMpass_IA20_STAND_ALONE_SURVIVAL.tsv")

# mmrf.pt.visit.with.seq$SPECTRUM_SEQ gives a list of potential samples to pull information from
# pts kept by JKeats in WES+WGS
keats.mut.full <- read_tsv("../../../MMRF_IA20/Somatic Mutation Files - SNV and INDEL_MMRF_CoMMpass_IA20_combined_vcfmerger2_All_Canonical_Variants.tsv", num_threads = 4)
keats.cna.full <- read_tsv("../../../MMRF_IA20/Copy Number Estimates_MMRF_CoMMpass_IA20_genome_gatk_cna.seg", num_threads = 4)
keats.sv.full <- read_tsv("../../../MMRF_IA20/Structural Event Files_MMRF_CoMMpass_IA20_genome_manta.tsv", num_threads = 4)
keats.tx.full <- read_tsv("../../../MMRF_IA20/SeqFISH Files_MMRF_CoMMpass_IA20_genome_tumor_only_mm_igtx_pairoscope.tsv", num_threads = 4)
keats.cna.thres.full <- read_tsv("../../../MMRF_IA20/SeqFISH Files_MMRF_CoMMpass_IA20_genome_gatk_cna_seqFISH.tsv", num_threads = 4)

length(unique(keats.mut.full$SAMPLE))
length(unique(keats.cna.full$SAMPLE))
length(unique(keats.sv.full$SAMPLE))

# pull data
rsp.maf <- read_tsv("../data/MMRF_TvsN_correct_CoveredNormals_NoSerial.aggregated.tsv", num_threads = 4)
jb.seg <- rbindlist(list(read_tsv("../data/commpass_500_1.seg"), read_tsv("../data/commpass_500_2_3.seg"))) 

# load SV workflow output
svcomm.part1.v2 <- read_tsv("../data/SV_filtered_commpass_hg38_500_1.v2.txt")
svcomm.part2.v2 <- read_tsv("../data/SV_filtered_commpass_hg38_500_2_3.v2.txt")

sv.workflow <- rbindlist(list(svcomm.part1.v2, svcomm.part2.v2))

mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb <- mmrf.pt.baseline.visit.with.seq %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos") %in% keats.mut.full$SAMPLE) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos") %in% keats.cna.full$SAMPLE) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos") %in% keats.sv.full$SAMPLE) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos") %in% keats.tx.full$SAMPLE) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos") %in% keats.cna.thres.full$SAMPLE) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM") %in% rsp.maf$sample) %>%
  filter(paste0(SPECTRUM_SEQ, "_BM_CD138pos_pair") %in% jb.seg$sample) %>%
  mutate(`Sample ID`=PUBLIC_ID)

mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb %>% write_tsv("../data/mmrf_final_annot.tsv")
mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb <- read_tsv("../data/mmrf_final_annot.tsv")

# export baseline svs -----------------------------------------------------

mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_pair_id <- paste0(mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$SPECTRUM_SEQ, "_BM_CD138pos_pair")

sv.workflow.baseline <- sv.workflow %>% 
  filter(individual %in% mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_pair_id)

length(unique(sv.workflow.baseline$individual))
#776
# SG Workflow likely no SV in few individuals (?)- will use IGH translocation forcecall and MMRF annotation to fill in blanks

sv.workflow %>% write_tsv("../data/20230614_MMRF_SV_filtered_results_hg38_n_982.tsv")
sv.workflow.baseline %>% write_tsv("../data/20230614_MMRF_SV_filtered_results_baseline_hg38_n_776.tsv")

# meta --------------------------------------------------------------------

# RSP Exome analysis ------------------------------------------------------
clinic.data <- read_tsv("../data/20230926_mmrf_clinical_data.tsv")

# now keep potential variants for matrix and QC (not exactly same list as in funcotator)
ns.variants <- c("Start_Codon_Del", "Stop_Codon_Del", "De_novo_Start_OutOfFrame",
                 "De_novo_Start_InFrame", "In_Frame_Ins", "Nonstop_Mutation",
                 "Start_Codon_SNP", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del",
                 "Splice_Site", "Nonsense_Mutation", "Missense_Mutation")
c.variants <- c("Silent", ns.variants)

all <- expand_grid(`Sample ID`=clinic.data$`Sample ID`, Hugo_Symbol=genes.of.interest)

clean.rsp.maf <- rsp.maf %>%
  filter(sample %in% paste0(paste0(clinic.data$SPECTRUM_SEQ, "_BM"))) %>%
  mutate(tumor_sample_barcode=sample)

# mut.counts <- clean.rsp.maf %>% group_by(sample, Variant_Classification) %>% count() %>% arrange(-n)

clean.rsp.maf %>% write_tsv("../data/20230503_export_mmrf_exome_maf_hg19_final_list.tsv")

#sanity check
genes.of.interest[which(genes.of.interest %nin% clean.rsp.maf$Hugo_Symbol)]
"H1-4" %in% clean.rsp.maf$Hugo_Symbol
"TENT5C" %in% clean.rsp.maf$Hugo_Symbol

# pull data
mut.bins <- clean.rsp.maf %>%
  filter(Hugo_Symbol %in% genes.of.interest) %>%
  filter(Variant_Classification %in% non.synonymous) %>%
  left_join(clinic.data, by=c("sample"="RSP_ID")) %>%
  right_join(all, by=c("Sample ID", "Hugo_Symbol")) %>%
  mutate(name=Hugo_Symbol) %>%
  group_by(`Sample ID`, name) %>% 
  summarise(value = prioritize.variant.classification(Variant_Classification[Variant_Classification!=""]),
            all_mutations = paste0(Variant_Classification[Variant_Classification!=""], collapse = ";", recycle0 = TRUE),
            n_mutations = length(Variant_Classification[Variant_Classification!="" & !is.na(Variant_Classification)])) %>% 
  mutate(value = case_when(is.na(value)|value=="NA"~"",TRUE~value)) %>% 
  mutate(all_mutations = case_when(is.na(all_mutations)|all_mutations=="NA"~"",TRUE~all_mutations)) %>%
  as.data.frame()

mut.bins %>% write_tsv("../data/fig1_tmp/mmrf_snv.tsv")

# CNV ---------------------------------------------------------------------

cnv.bins <- mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb %>%
  left_join(mmrf.pt, by=c("PUBLIC_ID")) %>%
  mutate(SEQ_SAMPLE_ID=paste0(SPECTRUM_SEQ, "_BM_CD138pos")) %>%
  left_join(keats.cna.thres.full, by=c("SEQ_SAMPLE_ID"="SAMPLE")) %>% 
  mutate(Del1p22.1 = case_when(SeqWGS_Cp_1p22_20percent==1~"loss", TRUE~""),
         Gain1q=case_when(SeqWGS_Cp_1q21_20percent==1~"gain", TRUE~""),
         Gain3q26.2=case_when(str_detect(SEQ_SAMPLE_ID, "2272|2307|2541|1781|2429|1839|1855|2476|1906|2539|2795|1787|1819|1842|2137|1550|1534|2365|2503|2378|2327|1722|1890|2005|1817|1955|2853|1453|2412")~"gain", TRUE~""),
         Del6q22.33=case_when(SeqWGS_Cp_6p22_20percent==1~"loss", TRUE~""),
         Del8p22=case_when(SeqWGS_Cp_8p22_20percent==1~"loss", TRUE~""),
         Del13q=case_when(as.numeric(SeqWGS_Cp_13q14_20percent+SeqWGS_Cp_13q34_20percent >=1 )==1~"loss", TRUE~""),
         Gain16p13.13=case_when(str_detect(SEQ_SAMPLE_ID, "1342|1957|2705|2801|1671|1650|1677|2279|2708|2595|2180|2232|2532|1633")~"gain", TRUE~""),
         Del14q32.32=case_when(SeqWGS_Cp_14q23_20percent==1~"loss", TRUE~""),
         Del17p=case_when(SeqWGS_Cp_17p13_20percent==1~"loss", TRUE~""),
         Del22q=case_when(SeqWGS_Cp_22q13_20percent==1~"loss", TRUE~"")) %>%
  select(`Sample ID`, Del1p22.1, Gain1q, Gain3q26.2, Del6q22.33, Del8p22, Del13q, Del14q32.32, Gain16p13.13, Del17p, Del22q) %>%
  pivot_longer(-`Sample ID`)

cnv.bins %>% write_tsv("../data/fig1_tmp/mmrf_cnv.tsv")

cnv.bins %>%
  ggplot(aes(`Sample ID`, name, fill=value)) +
  geom_raster()


# JB genomes analysis -----------------------------------------------------
# wolf gneomes hg38, liftOvered hg19
# not exome
# MAY 3: check again that it is working fine
# Not retested. Trust my code.
# commpass_maf_sample_hg19.maf minimal set of variables

# mmrf.genomes <- read_tsv("../maf-for-merge-ctc-su2c/commpass_maf_sample_hg19.maf", num_threads = 4)
mmrf.genomes <- read_tsv("../maf-for-merge-ctc-su2c/20230217_mmrf_commpass_baseline_hg19.maf", num_threads = 4)

mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb <- read_tsv("../data/mmrf_final_annot.tsv")
mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_pair_id <- 
  paste0(mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$SPECTRUM_SEQ, "_BM_CD138pos_pair_tumor")
mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_ctf_id <- 
  paste0(mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$SPECTRUM_SEQ, "_BM_CD138pos")

# mmrf.genomes$Tumor_Sample_Barcode
# in the form
# MMRF_2041_1_BM_CD138pos_pair_tumor

head(mmrf.genomes$Tumor_Sample_Barcode)
head(mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_pair_id)

# subtract only first sample of BM
# 1726/1715 have clear indel problem in lib prep

mmrf.genomes.baseline <- mmrf.genomes %>% filter(Tumor_Sample_Barcode %in% mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_pair_id)
mmrf.genomes.baseline <- mmrf.genomes.baseline %>% filter(Tumor_Sample_Barcode != "MMRF_1726_1_BM_CD138pos_pair_tumor")
mmrf.genomes.baseline <- mmrf.genomes.baseline %>% filter(Tumor_Sample_Barcode != "MMRF_1715_1_BM_CD138pos_pair_tumor")

gc()
# 776 patients final
#
# mmrf.genomes.baseline %>% 
#   group_by(Tumor_Sample_Barcode, Variant_Type) %>%
#   count() %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   summarise(fraction = n[Variant_Type=="INS"]/n[Variant_Type=="SNP"]) %>%
#   View
# mmrf.genomes.baseline %>% 
#   group_by(Tumor_Sample_Barcode, Variant_Classification) %>%
#   count() %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   summarise(fraction = n[Variant_Classification=="Silent"]/n[Variant_Classification=="Missense_Mutation"]) %>%
#   View
# table(mmrf.genomes.baseline$Tumor_Sample_Barcode)

# filter(SPECTRUM_SEQ %in% mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$SPECTRUM_SEQ)

cols.to.keep <- c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", 
  "Start_position", "End_position", "Strand", "Variant_Classification", 
  "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
  "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Genome_Change", "Protein_Change",
  "t_alt_count", "t_ref_count", "ccf_hat", "ccf_CI_low", "ccf_CI_high", 
  "ref_context")
mmrf.genomes.baseline.light <- mmrf.genomes.baseline %>% select(all_of(cols.to.keep))
mmrf.genomes.light <- mmrf.genomes %>% select(all_of(cols.to.keep))

mmrf.genomes %>% write_tsv("../data/20230505_mmrf_commpass_hg19.maf.gz", num_threads = 4)
mmrf.genomes.light %>% write_tsv("../data/20230505_mmrf_commpass_hg19_light.maf.gz", num_threads = 4)
mmrf.genomes.baseline %>% write_tsv("../data/20230505_mmrf_commpass_baseline_hg19.maf.gz", num_threads = 4)
mmrf.genomes.baseline.light %>% write_tsv("../data/20230505_mmrf_commpass_baseline_hg19_light.maf.gz", num_threads = 4)
gc()

# TX analysis -------------------------------------------------------------

tx.bins <- mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb %>%
  left_join(mmrf.pt, by=c("PUBLIC_ID")) %>%
  mutate(SEQ_SAMPLE_ID=paste0(SPECTRUM_SEQ, "_BM_CD138pos")) %>%
  mutate(RSP_ID=paste0(SPECTRUM_SEQ, "_BM")) %>%
  left_join(keats.tx.full, by=c("SEQ_SAMPLE_ID"="SAMPLE"))

# trick to multiple CALL (0 or 1) with IGSOURCE (1, 2, 3) to eliminate default source on no-call
rename.columns.tx.keats <- tx.bins %>%
  select(SEQ_SAMPLE_ID, any_of( c(ends_with("CALL", ignore.case=FALSE), ends_with("IGSOURCE")) )) %>%
  pivot_longer(cols = -c(SEQ_SAMPLE_ID), names_to = c(".value", "gene"), names_pattern = "(.*)_(.*)") %>%
  group_by(SEQ_SAMPLE_ID) %>%
  summarise(across(-gene, prod)) %>%
  group_by(SEQ_SAMPLE_ID) %>%
  mutate(across(everything(), ~ case_when(.x=="1" ~ "IGH", .x=="2" ~ "IGK", .x=="3" ~ "IGL", TRUE ~ NA))) %>%
  rename("MMSET"="NSD2")

keats.tx.ig.calls <- rename.columns.tx.keats |>
  pivot_longer(cols=c(-SEQ_SAMPLE_ID), names_to = "onco", values_to = "ig") |>
  filter(!is.na(onco)&!is.na(ig)) |>
  mutate(igonco=paste(ig, onco)) |>
  select(-onco, -ig) |>
  mutate(keep=1) |> 
  group_by(igonco) |>
  mutate(count=n()) |>
  arrange(-count) |> 
  pivot_wider(id_cols = c(SEQ_SAMPLE_ID), names_from = igonco, values_from = keep)
  # pivot_wider(id_cols = c(SEQ_SAMPLE_ID, ), names_from = "onco", names_glue = "{.ig} {onco}")
# 

# TODO 
#
# FIX TX BINS
#
# tx.bins %>% write_tsv("../data/fig1_tmp/mmrf_tx.tsv")
# 
# tx.bins %>%
#   ggplot(aes(`Sample ID`, name, fill=value)) +
#   geom_raster()


# SV_WORKFLOW -------------------------------------------------------------

ctf.ig.coordinates <- read_tsv("~/catchthefish/ref/GRCh38-ig-regions-CoMMPASS-adjusted.bed", col_names = c("chr", "start", "end", "ig"))
ctf.onco.coordinates <- read_tsv("~/catchthefish/ref/GRCh38-oncogenes-regions-CoMMPASS-adjusted.bed", col_names = c("chr", "start", "end", "onco"))

target.candidates.sv.workflow.drivers <- sv.workflow.baseline |>
  mutate(chr1=as.character(chr1), chr2=as.character(chr2)) |>
  filter(class=="inter_chr") |>
  left_join(ctf.ig.coordinates, by=c("chr1"="chrN")) |>
  mutate(ig_final=case_when(
    pos1>=start&pos1<=end~ig, 
    TRUE~NA)) |>
  select(-chr, -start, -end, -ig) |>
  left_join(ctf.ig.coordinates, by=c("chr2"="chrN")) |>
  mutate(ig_final=case_when(
    (!is.na(ig_final)) & ( pos2>=start&pos2<=end) ~ NA, #don't allow IG-IG reports 
    (is.na(ig_final)) & ( pos2>=start&pos2<=end) ~ ig, 
    TRUE~ig_final)) |>
  select(-chr, -start, -end, -ig) |> 
  left_join(ctf.onco.coordinates, by=c("chr1"="chrN"), relationship = "many-to-many") |> # expected when potentially many targets per chromosome, next row will keep max 1 candidate
  mutate(onco_final=case_when(
    pos1>=start&pos1<=end~onco, 
    TRUE~NA)) |>
  select(-chr, -start, -end, -onco) |>
    left_join(ctf.onco.coordinates, by=c("chr2"="chrN"), relationship = "many-to-many") |>
    mutate(onco_final=case_when(
      (!is.na(onco_final)) & ( pos2>=start&pos2<=end) ~ paste(onco_final, onco), # allow ONCO-ONCO reports 
      (is.na(onco_final)) & ( pos2>=start&pos2<=end) ~ onco, 
      TRUE~onco_final)) |>
  select(-chr, -start, -end, -onco) |>
  filter(!is.na(onco_final) & (!is.na(ig_final) | str_detect(onco_final, " "))) |>
  mutate(ig_onco=case_when(is.na(ig_final) ~ onco_final,
                           TRUE ~ paste(ig_final, onco_final))) |>
  relocate(ig_final, onco_final, ig_onco)

sv.workflow.ig.calls <- target.candidates.sv.workflow.drivers |>
  select(ig_onco, individual) |> 
  distinct() |>
  group_by(ig_onco) |>
  mutate(count=n()) |>
  arrange(-count) |>
  mutate(keep=1) |>
  pivot_wider(id_cols = individual, names_from = ig_onco, values_from = keep) |>
  mutate(individual=str_remove(individual, "_pair$"))

# CTF ---------------------------------------------------------------------

ctf <- read_tsv("../../WGS_JB/MMRF_WGS_1933_samples.aggregated.tsv")

manual.inspection.ids <- c("MMRF_1269_1_BM_CD138pos", "MMRF_1439_1_BM_CD138pos", "MMRF_1647_1_BM_CD138pos", 
                           "MMRF_1817_1_BM_CD138pos", "MMRF_2014_1_BM_CD138pos", "MMRF_2017_1_BM_CD138pos", 
                           "MMRF_2308_1_BM_CD138pos", "MMRF_2472_1_BM_CD138pos", "MMRF_1683_1_BM_CD138pos", 
                           "MMRF_1361_1_BM_CD138pos", "MMRF_1391_1_BM_CD138pos", "MMRF_1572_1_BM_CD138pos", 
                           "MMRF_1823_1_BM_CD138pos", "MMRF_2477_1_BM_CD138pos", "MMRF_1432_1_BM_CD138pos", 
                           "MMRF_2326_1_BM_CD138pos", "MMRF_2153_1_BM_CD138pos", "MMRF_1403_1_BM_CD138pos")

ctf %>% filter(startsWith(ID, manual.inspection.ids)) %>% View()

ctf.ig.calls <- ctf |> 
  filter(INTERVAL_LENGTH>50 & IG_LENGTH>50) |> 
  filter(maxMAPQ.ONCO>=59) |> 
  filter(N_Molecules>=3) |>
  mutate(individual=str_extract(ID, "^MMRF.*CD138pos")) |>
  filter(individual %in% mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb$expected_ctf_id) |>
  relocate(individual) |>
  select(-ID) |>
  group_by(individual) |>
  mutate(IGONCO=paste(IG, ONCO)) |>
  select(individual, IGONCO) |>
  distinct() |>
  group_by(IGONCO) |>
  mutate(count=n()) |>
  arrange(-count) |>
  mutate(keep=1) |>
  pivot_wider(id_cols = individual, names_from = IGONCO, values_from = keep)


# Aggregate tx ------------------------------------------------------------

# aggreg from Keats, SV_WORKFLOW, and catchthefish

ctf.ig.calls %>%
  full_join(sv.workflow.ig.calls, by=c("individual"="individual"), suffix = c(" CTF", " SVW")) %>%
  full_join(keats.tx.ig.calls %>% rename_all(~paste0(.x, " KTS")), by=c("individual"="SEQ_SAMPLE_ID KTS"), suffix = c("", " KTS")) %>%
  select(sort(names(.))) %>%
  relocate(individual) %>%
  clipr::write_clip()

final.review.tx <- read_xlsx("../data/MMRF_ANNOT_TX_SV_CTF_KTS.xlsx")

sort(table(final.review.tx$`Primary Tx Final Call`), decreasing = TRUE)


# META --------------------------------------------------------------------

# putting this down since tx come from meta-analysis of ctf wolf and keats

clinic.data <- mmrf.pt.baseline.with.wxs.and.wgs.and.rsp.and.jb %>%
  left_join(final.review.tx, by=c("expected_ctf_id"="individual")) %>%
  left_join(mmrf.pt, by=c("PUBLIC_ID")) %>%
  mutate(SEQ_SAMPLE_ID=paste0(SPECTRUM_SEQ, "_BM_CD138pos")) %>%
  mutate(RSP_ID=paste0(SPECTRUM_SEQ, "_BM")) %>%
  left_join(keats.cna.thres.full, by=c("SEQ_SAMPLE_ID"="SAMPLE")) %>%
  mutate(Gender=case_when(D_PT_gender==1~"M", 
                          D_PT_gender==2~"F",
                          TRUE~"O"),
         Stage="NDMM",
         Age=D_PT_age,
         HRDTx=case_when(!is.na(`Primary Tx Final Call`)~`Primary Tx Final Call`,
                         SeqWGS_Cp_Hyperdiploid_Call==1 ~ "HRD",
                         TRUE ~ "UNKNOWN"),
         Cohort="MMRF",
         Assay="W(G+X)S") %>%
  relocate(SPECTRUM_SEQ, Gender, Stage, Age, HRDTx, Cohort, Assay)

clinic.data %>% write_tsv("../data/20230926_mmrf_clinical_data.tsv")
clinic.data <- read_tsv("../data/20230926_mmrf_clinical_data.tsv")

clinic.bins <- clinic.data %>%
  mutate(Age=as.character(Age)) %>%
  select(`Sample ID`, Gender, Age, Stage, HRDTx, Cohort, Assay) %>%
  pivot_longer(-`Sample ID`)

clinic.bins %>% write_tsv("../data/fig1_tmp/mmrf_meta.tsv")

primary.tx <- clinic.data %>%
  select(`Sample ID`, `Primary Tx Final Call`) %>%
  mutate(value=case_when(!is.na(`Primary Tx Final Call`) | `Primary Tx Final Call`=="" ~"tx", TRUE~"")) %>%
  pivot_wider(names_from = 2) %>%
  select(-`NA`) %>%
  pivot_longer(-`Sample ID`)
  
secondary.tx <- clinic.data %>%
  select(`Sample ID`, `Secondary Tx Final Call`) %>%
  mutate(name="t(MYC)",
         value=case_when(str_detect(`Secondary Tx Final Call`, "MYC") ~ "tx", TRUE~"")) %>%
  select(-`Secondary Tx Final Call`)

tx.bins <- rbind(primary.tx, secondary.tx)

tx.bins %>% write_tsv("../data/fig1_tmp/mmrf_tx.tsv")

#   
# tx.bins <- clinic.data %>%
#   select(`Sample ID`, `Primary Tx Final Call`, `Secondary Tx Final Call`) %>%
#   mutate(`Secondary Tx Final Call`=case_when( str_detect(`Secondary Tx Final Call`, "MYC") ~ "t(MYC)", TRUE~""),
#          `Primary Tx Final Call`=case_when(!is.na(`Primary Tx Final Call`)~`Primary Tx Final Call`, TRUE~"")) %>%
#   pivot_longer(-`Sample ID`)

# 
# tx.bins |> head()
#   pivot_wider(id_cols = c(`Sample ID`, `Primary Tx Final Call`), names_from = "value") |> 
#   pivot_wider(id_cols = c(`Sample ID`, ), names_from = "value") |> View()
#   select(-HRD, -UNKNOWN) |> View()
#   mutate(across(starts_with("t("), ~case_when(!is.na(.x) ~ "tx", TRUE ~ ""))) |>
#   pivot_longer(-`Sample ID`) |> View()
#   write_tsv("../data/fig1_tmp/mmrf_tx.tsv")

clinic.bins %>% write_tsv("../data/fig1_tmp/mmrf_meta.tsv")

clinic.bins %>%
  ggplot(aes(`Sample ID`, name, fill=value)) +
  geom_raster()
