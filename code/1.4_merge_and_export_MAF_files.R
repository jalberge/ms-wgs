# This merged maf will be used to run MutSig2CV and annotate patients
# ( note from July 10 2024, Some RedCap IDs have been hidden here as well and cannot be made public  - to be changed with dbGAP ID)
# !! see here that indels with coverage (REF+ALT) < 8 have been removed.
# This is to homogenize with MuTect rule that ALT+REF should be at least 8. 
# This REF+ALT rule is missing in wolF workflow for Strelka indels

library(tidyverse)
library(readxl)

JCO.SMM.MAF <- "../data/20230509_bustoros_exome_hg19.maf"
MMRF.WXS.RSP.MAF <- "../data/20230503_export_mmrf_exome_maf_hg19_final_list.tsv"
SU2C.PT.SPECIFIC.MAF <- "../data/20230503_su2c_oben_wtc_ctc_wm_prm_n_182_patient_level_hg19.tsv.gz"
SU2C.FULL.MAF <- "../data/20230503_su2c_oben_wtc_ctc_wm_prm_n_228_sample_level_hg19.tsv.gz"
MMRF.GENOMES.MAF <- "../data/20230505_mmrf_commpass_baseline_hg19.maf.gz"

# GENOME AND EXOME --------------------------------------------------------

maf.full.hg19 <- read_tsv(SU2C.FULL.MAF, num_threads = 4)
maf.full.hg19.patient.specific <- read_tsv(SU2C.PT.SPECIFIC.MAF, num_threads = 4)
mark.jco <- read_tsv(JCO.SMM.MAF, num_threads = 4)
rsp.mmrf <- read_tsv(MMRF.WXS.RSP.MAF, num_threads = 4)

mark.jco$NCBI_Build <- "hg19"
mark.jco$Tumor_Sample_Barcode=mark.jco$sample

rsp.mmrf$NCBI_Build <- "hg19"
rsp.mmrf$Tumor_Sample_Barcode=rsp.mmrf$sample


# Exomes/Genomes ----------------------------------------------------------

mark.jco <- mark.jco[, colnames(mark.jco) %in% colnames(maf.full.hg19)]
rsp.mmrf <- rsp.mmrf[, colnames(rsp.mmrf) %in% colnames(maf.full.hg19)]

maf.full.hg19.for.mmrf.jco <- maf.full.hg19.patient.specific[, colnames(maf.full.hg19.patient.specific) %in% colnames(mark.jco)]

su2c.mark.jco.rsp.mmrf <- rbindlist(list(maf.full.hg19.for.mmrf.jco, mark.jco, rsp.mmrf), fill = TRUE)

head(sort(unique(su2c.mark.jco.rsp.mmrf$Tumor_Sample_Barcode)))
tail(sort(unique(su2c.mark.jco.rsp.mmrf$Tumor_Sample_Barcode)))

su2c.mark.jco.rsp.mmrf <- su2c.mark.jco.rsp.mmrf %>%
  filter(!(Tumor_Sample_Barcode %in% c("xxxx", "xxxx"))) %>%
  filter(Chromosome %in% c(1:22, "X", "Y"))

write_tsv( su2c.mark.jco.rsp.mmrf, "../data/20230509.su2c.mark.jco.rsp.mmrf.hg19.nowm.nomgip.maf")

su2c.mark.jco.rsp.mmrf <- read_tsv("../data/20230509.su2c.mark.jco.rsp.mmrf.hg19.nowm.nomgip.maf", num_threads = 4)

# GENOMES -----------------------------------------------------------------

maf.full.hg19.patient.specific <- read_tsv(SU2C.PT.SPECIFIC.MAF, num_threads = 4)
maf.full.hg19.patient.specific <- maf.full.hg19.patient.specific %>%
  mutate(Participant_ID=case_when(is.na(Participant_ID)~Tumor_Sample_Barcode, TRUE~Participant_ID))

maf.full.hg19 <- read_tsv(SU2C.FULL.MAF, num_threads = 4)
maf.full.hg19$Participant_ID <- maf.full.hg19$Patient

maf.hg19.commpass.baseline <- read_tsv(MMRF.GENOMES.MAF, num_threads = 4)
maf.hg19.commpass.baseline$Participant_ID <- str_extract(maf.hg19.commpass.baseline$Tumor_Sample_Barcode, "MMRF_[0-9]{4}")

maf.full.hg19.patient.specific <- maf.full.hg19.patient.specific[, colnames(maf.full.hg19.patient.specific) %in% colnames(maf.hg19.commpass.baseline)]
maf.full.hg19 <- maf.full.hg19[, colnames(maf.full.hg19) %in% colnames(maf.hg19.commpass.baseline)]

tmp.mmrf.genomes <- maf.hg19.commpass.baseline[, colnames(maf.hg19.commpass.baseline) %in% colnames(maf.full.hg19.patient.specific)]
tmp.mmrf.genomes.2 <- maf.hg19.commpass.baseline[, colnames(maf.hg19.commpass.baseline) %in% colnames(maf.full.hg19)]

su2c.mmrf.genomes <- data.table::rbindlist(list(maf.full.hg19.patient.specific, tmp.mmrf.genomes), use.names = TRUE)
su2c.mmrf.genomes.2 <- data.table::rbindlist(list(maf.full.hg19, tmp.mmrf.genomes.2), use.names = TRUE)

# clean dataset -----------------------------------------------------------

su2c.mmrf.genomes %>%
  group_by(Variant_Type) %>%
  summarise(FILTER=sum( n_alt_count+n_ref_count < 8), N=n())

su2c.mmrf.genomes %>%
  group_by(Variant_Type) %>%
  summarise(FILTER=sum( t_alt_count+t_ref_count < 6), N=n())

su2c.mmrf.genomes.clean <- su2c.mmrf.genomes %>%
  filter(Chromosome %in% c(1:22, "X", "Y")) %>% #lifted over realignment
  filter(n_alt_count+n_ref_count >= 8) # MuTect rules not forwarded to Strelka in the wgs workflow

su2c.mmrf.genomes.clean.sample.level <- su2c.mmrf.genomes.2 %>%
  filter(Chromosome %in% c(1:22, "X", "Y")) %>% #lifted over realignment
  filter(n_alt_count+n_ref_count >= 8) # MuTect rules not forwarded to Strelka in the wgs workflow

# export for signature analyzer / hdp / dig
su2c.mmrf.genomes.clean %>% write_tsv("../data/20230614_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_patient_level.tsv.gz", num_threads = 4)
su2c.mmrf.genomes.clean.sample.level %>% write_tsv("../data/20230614_su2c_mmrf_genomes_hg19_1-22-XY_cleanIndels_light_sample_level.tsv.gz", num_threads = 4)


