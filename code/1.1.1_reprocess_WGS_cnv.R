# CNV Reprocessing of WGS hg38 and hg19 (lifterOver) for Andrea to run GISTIC and League Model

# INIT --------------------------------------------------------------------
source("code/0_annotate_samples.R")

# Load CNV ----------------------------------------------------------------

# ABSOLUTE output files
# seg files kept integer notation for contigs even with hg38 ref
seg.su2c <- read_tsv("../data/SU2C_hg38_20230130.seg")
seg.ctc <- read_tsv("../data/CTC_hg19_20230130.seg")
seg.mb <- read_tsv("../data/MB_hg38_20230130.seg")
seg.obenwtc <- read_tsv("../data/OBENWTC_hg38_20230130.seg")

# seg.su2c$REF="hg38"
# seg.mb$REF="hg38"
# seg.obenwtc$REF="hg38"
# seg.ctc$REF="hg19"

# RUN GISTIC 2 SU2C
# RUN GISTIC 2 COMMpass
# VALIDATE CTC hg19 (MinimuMM-seq paper)
# ANNOTATE

# for any hg19 project, use liftOverWrapperSEG in maf-for-merge-ctc-su2c

seg.ctc.hg38 <- read_tsv("../data/CTC_hg19_20230130_liftOveredhg38.seg")
#because of liftOver, needs to remove ALT and chr prefix
seg.ctc.hg38 <- seg.ctc.hg38 %>%
  filter(Chromosome %in% paste0("chr", 1:22)) %>%
  mutate(Chromosome=as.numeric(str_remove(Chromosome, "chr"))) %>%
  arrange(sample, Chromosome, Start.bp)

# remove overlapping fragments (undecided CN?)
seg.ctc.hg38 <- seg.ctc.hg38 %>%
  group_by(sample, Chromosome) %>%
  filter( ! (End.bp > lead(End.bp) & !is.na(lead(End.bp))) )
  
seg.ctc.hg38 %>% group_by(sample, Chromosome) %>% filter(End.bp==lead(End.bp))
# none

# remaining: a couple of sliced pericentromeric regions of chr1
seg.ctc.hg38 <- seg.ctc.hg38 %>%
  group_by(sample, Chromosome) %>% 
  mutate(End.bp = ifelse( End.bp > lead(Start.bp) & !is.na(lead(End.bp)), lead(Start.bp) , End.bp))

# Add log2 space scale
all.seg.files <- rbindlist(list(seg.su2c, seg.mb, seg.obenwtc, seg.ctc.hg38)) %>%
  mutate(Seg.CN = log2( rescaled_total_cn ) -1)

all.segs <- sort(unique(all.seg.files$sample))

all.segs[all.segs %nin% clinical.and.terra.ref.pairs.samples$Reference_Pair]
  
# FIRST EXPORT ALL --------------------------------------------------------

all.seg.files %>%
  select(sample, Chromosome, Start.bp, End.bp, n_probes, Seg.CN) %>%
  write_tsv("../data/hg38_all_wgs_for_gistic.seg")

# Then export 1 per patient -----------------------------------------------

all.seg.files %>%
  filter(sample %in% clinical.and.terra.ref.pairs.samples$Reference_Pair) %>%
  select(sample, Chromosome, Start.bp, End.bp, n_probes, Seg.CN) %>%
  write_tsv("../data/hg38_reference_pair_wgs_for_gistic.seg")

# Then export only MGUS/SMM -----------------------------------------------

all.seg.files %>%
  filter(sample %in% clinical.and.terra.ref.pairs.samples$Reference_Pair[clinical.and.terra.ref.pairs.samples$Disease_Status %in% c("HRSMM", "IRSMM", "LRSMM", "MGUS", "SMM")]) %>%
  select(sample, Chromosome, Start.bp, End.bp, n_probes, Seg.CN) %>%
  write_tsv("../data/hg38_reference_pair_mgus_smm_only_wgs_for_gistic.seg")

