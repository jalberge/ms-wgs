# README
# THIS SCRIPT TO LOAD LYMPH B / T SINGLE COLONY SEQUENCING FROM SANGER
# https://github.com/machadoheather/lymphocyte_somatic_mutation/blob/main/data/som2_indel_SNV_allDF_v2.txt.gz
# https://www.nature.com/articles/s41586-022-05072-7

# load sanger mut ---------------------------------------------------------

library(tidyverse)
library(stringr)

setwd("code/")
source("0_annotate_samples.R")

HOMEMACHADO <- "~/Dropbox (Partners HealthCare)/lymphocyte_somatic_mutation/"

hmach <- read_tsv(paste0(HOMEMACHADO, "data/som2_indel_SNV_allDF_v2.txt.gz"))
colonyinfo <- read_tsv(paste0(HOMEMACHADO, "data/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt"))
mutcount <- read_tsv(paste0(HOMEMACHADO, "data/mutcounts_matrix_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt"))

# ggplot(hmach, aes(HQ.alt.freq, color=sample)) + geom_density()

b.hmach <- hmach %>%
  left_join(colonyinfo, by=c("sample"="colony")) %>%
  filter(CellType %in% c("B", "B Memory", "B Naive"))
  
nrow(b.hmach)
head(b.hmach)
# first subset normal b cells
# subfirst match patient sample to phenotype
# second vcf to maf
# subsecond need to recode position start / end for INDEL in particular. Also check +1/-1

simple.b.hmach <- b.hmach %>%
  filter(nchar(ref)==1 & nchar(mut)==1)

complex.b.hmach <- b.hmach %>% filter( nchar(b.hmach$ref) >= 2 & nchar(b.hmach$mut) >= 2) 
simple.b.hmach <- b.hmach %>% filter( nchar(b.hmach$ref) < 2 | nchar(b.hmach$mut) < 2) 

simple.b.hmach %>% filter(nchar(ref)>=2) %>% View
simple.b.hmach %>% filter(nchar(mut)>=2) %>% View

simple.b.hmach <- simple.b.hmach %>%
  mutate(vt=case_when( nchar(ref) == 1 & nchar(mut) == 1 ~ "SNP",
                       nchar(ref) >= 2 & nchar(mut) == 1 ~ "DEL",
                       nchar(ref) == 1 & nchar(mut) >= 2 ~ "INS",
                       TRUE ~ "ERROR") )


# compare with our maf ----------------------------------------------------


simple.b.hmach.maf <- simple.b.hmach %>%
  mutate(Hugo_Symbol="tbd",
         Entrez_Gene_Id="tbd",
         Center="Sanger",
         NCBI_Build="hg19",
         Chromosome=chr,
         Start_position=case_when(vt=="DEL"~pos+1, TRUE~pos),
         End_position=case_when(vt=="DEL"~pos+nchar(ref)-1, vt=="INS"~pos+1,TRUE~pos),
         Strand=".",
         Variant_Classification=".",
         Variant_Type=vt,
         Reference_Allele=case_when(vt=="DEL"~str_sub(ref, 2, -1), vt=="INS"~"-", TRUE~ref),
         Tumor_Seq_Allele1=Reference_Allele,
         Tumor_Seq_Allele2=case_when(vt=="DEL"~"-", vt=="INS"~str_sub(mut, 2, -1), TRUE~mut),
         Tumor_Sample_Barcode=sample,
         Matched_Norm_Sample_Barcode="UNKNOWN",
         .before = ID)

simple.b.hmach.maf %>% write_tsv("../data/20230228_try_vcf2maf_conv_normal_B_colonies_Machado.maf")
