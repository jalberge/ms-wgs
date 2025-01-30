
# Packages ----------------------------------------------------------------

library(tidyverse)
library(data.table)
library(rstatix)
library(kableExtra)
library(RColorBrewer)
# library(maftools)
library(cowplot)
library(readxl)
library(ggpubr)
library(scales)
library(circlize)
library(colorspace)

# Variables ---------------------------------------------------------------

# these breaks are mostly used for pseudo_log scales in ggplot
int.breaks <- c(0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000)

breaks <- c(1e-10, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 0.001, 0.01, 
            0.1, 1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 
            1e+10)
minor_breaks <- c(1e-10, 2e-10, 3e-10, 4e-10, 5e-10, 6e-10, 7e-10, 8e-10, 9e-10, 
                  1e-09, 2e-09, 3e-09, 4e-09, 5e-09, 6e-09, 7e-09, 8e-09, 9e-09, 
                  1e-08, 2e-08, 3e-08, 4e-08, 5e-08, 6e-08, 7e-08, 8e-08, 9e-08, 
                  1e-07, 2e-07, 3e-07, 4e-07, 5e-07, 6e-07, 7e-07, 8e-07, 9e-07, 
                  1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06, 7e-06, 8e-06, 9e-06, 
                  1e-05, 2e-05, 3e-05, 4e-05, 5e-05, 6e-05, 7e-05, 8e-05, 9e-05, 
                  1e-04, 2e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 8e-04, 9e-04, 
                  0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 
                  0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 
                  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                  10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 
                  600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 
                  8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 
                  80000, 90000, 1e+05, 2e+05, 3e+05, 4e+05, 5e+05, 6e+05, 7e+05, 
                  8e+05, 9e+05, 1e+06, 2e+06, 3e+06, 4e+06, 5e+06, 6e+06, 7e+06, 
                  8e+06, 9e+06, 1e+07, 2e+07, 3e+07, 4e+07, 5e+07, 6e+07, 7e+07, 
                  8e+07, 9e+07, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08, 6e+08, 7e+08, 
                  8e+08, 9e+08, 1e+09, 2e+09, 3e+09, 4e+09, 5e+09, 6e+09, 7e+09, 
                  8e+09, 9e+09, 1e+10, 2e+10, 3e+10, 4e+10, 5e+10, 6e+10, 7e+10, 
                  8e+10, 9e+10)

contigs <- factor(c(1:22, "X", "Y"), levels = c(1:22, "X", "Y"))

tau.max <- 4
tau.min <- 0


genes.of.interest.table <- read_tsv("../annot/20240619_mutsig2cv_more_than_1pct.txt", col_names = TRUE)
genes.of.interest.rare.table <- read_tsv("../annot/20240619_mutsig2cv_less_than_1pct.txt", col_names = TRUE)
genes.of.interest <- genes.of.interest.table %>% pull(1)

# grey for BMPCs, light blue for CMMCs
paired.pal=c(CMMC="#4EB3D3", BMPC="#B3B3B3")
paired.pals= c(CMMCs="#4EB3D3", BMPCs="#B3B3B3")

non.synonymous <- c("COULD_NOT_DETERMINE", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", "De_novo_Start_OutOfFrame",
                    "De_novo_Start_InFrame", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins", "In_Frame_Ins",
                    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",  "Splice_Site", "Start_Codon_SNP", "START_CODON_SNP", "Translation_Start_Site", "Start_Codon_Del", "Stop_Codon_Del")
coding <- c(non.synonymous, "Silent")

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# Functions ---------------------------------------------------------------

`%nin%` = Negate(`%in%`)

`%+%` <- paste0 # see discussions on overwriting + functions for chr https://stackoverflow.com/questions/4730551/making-a-string-concatenation-operator-in-r
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 1 is 4 copies (threshold)
# 1.3 is 5 copies approx
col_fun_tau = colorRamp2(c(0, 1.85, 2.15, 4), c("#2166AC", "#F7F7F7","#F7F7F7", "#B2182B"))
col_fun_l2r = colorRamp2(c(-.6, 0, .6), c("blue", "white", "red"))

# col_call = function(x) { ifelse(x=='+', '#FF0000', ifelse(x=='-', '#0000FF', "#FFFFFF"))}

# in cases where we have to annotate variants, pick first one by order or importance
prioritize.variant.classification <- function(x) {
  snv.order <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                 "Nonstop_Mutation", "Splice_Site", "START_CODON_SNP",
                 "Start_Codon_Del", "Stop_Codon_Del", 
                 "Start_Codon_SNP", "In_Frame_Del",  "In_Frame_Ins", "Missense_Mutation", 
                 "De_novo_Start_OutOfFrame", "De_novo_Start_InFrame",
                 "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", 
                 "Translation_Start_Site", "COULD_NOT_DETERMINE", "Silent")
  if(!any(is.na(x)))  as.character(head(sort(factor(x, levels=snv.order)), n=1))
  else NA
}

read_maf = function(maf, id=NULL) {
  dt <- data.table::fread(file = maf, sep = "\t", stringsAsFactors = FALSE, 
                    verbose = FALSE, data.table = TRUE, showProgress = TRUE, 
                    header = TRUE, fill = TRUE, skip = "Hugo_Symbol", 
                    quote = "")
  # changing tumor id is useful when mutation validator
  # (or other tools) join mutation without editing all the 
  # tumor_sample_barcode of origin
  if(!is.null(id)) dt$Tumor_Sample_Barcode=Tumor_Sample_UUID=id
  dt |> mutate_at("Entrez_Gene_Id", as.character)
}

load.m1fc <- function(power.wig, call_stats) {
  # takes MuTect1 force calling input files.
  # todo add coverage file
  pwig <- read_delim(power.wig, skip=2, col_names=c("Power"), delim="\t", comment="fixedStep")
  cs <-  read_delim(call_stats, delim="\t", comment="##")
  df <- cbind(pwig, cs)
  df
}

calculate_log_likelihood <- function(depth, alts, eps, f){
  a = (depth-alts) * log10( f*eps + (1-f)*(1-eps))
  b = alts * log10( f * (1-eps) + (1-f)*eps)
  a + b
}

calculate_t_lod <- function(depth, alts, eps=0.001, contam=0){
  f = alts/depth
  tlod = calculate_log_likelihood(depth, alts, eps, f) - calculate_log_likelihood(depth, alts, eps, pmin(f, contam))
  tlod
}

calculate_power <- function(depth, eps, lodThreshold=6.3, delta, contam=0) {
  
  if (depth == 0) return(0)
  
  p_alt_given_e_delta = delta * (1-eps) + (1-delta) * eps
  
  probs = dbinom(x = 0:depth, size = depth, prob = p_alt_given_e_delta)
  
  lods = calculate_t_lod(depth = depth, alts = 0:depth, eps = eps, contam = 0)
  
  k = max(0, min(which(lods>=lodThreshold))) # k is N-alt+1
  
  if(k > 0) {
    x = 1 - ( lodThreshold - lods[k-1] )/( lods[k] - lods[k-1])
    p0 = x * probs[k-1]
  } else {
    p0 = 0
  }
  
  power = p0 + sum(probs[k:length(probs)])
  power
}

calculate_tumor_power <- function(depth, eps, tlod, rho, psi, CCF, contam) {
  # adapted from MuTect1 java github code
  # TODO add tumor ploidy vs normal ploidy
  f = rho / psi * CCF
  calculate_power(depth = depth, eps = eps, lodThreshold = tlod, delta = f, contam = contam)
}

read_mixcr_clones <- function(x){
  read_tsv(x) |> filter(cloneId!="cloneId") |> mutate(isotype=str_extract(string=allVHitsWithScore, pattern="IG[HKL]"))
}

sensitivity <- function(real, measured) {
  real.pos = sum(real==1)
  true.measured.pos = sum(real==1 & measured==1)
  sensitivity = true.measured.pos / real.pos
  sensitivity
}

specificity <- function(real, measured) {
  real.neg = sum(real==0)
  true.measured.neg = sum(real==0 & measured==0)
  specificity = true.measured.neg / real.neg
  specificity
}

extract_first_vdj <- function(x) str_split(x, pattern = "[(]", )[[1]][1]

merge_maf_pair = function(maf1, maf2, suffix1, suffix2, id=NULL) {
  dt <- inner_join(
    maf1,
    maf2,
    by = c(
      "Hugo_Symbol",
      "Entrez_Gene_Id",
      "Center",
      "NCBI_Build",
      "Chromosome",
      "Start_position",
      "End_position",
      "Strand",
      "Variant_Classification", 
      "Variant_Type",
      "Reference_Allele"),
    suffix=c(suffix1, suffix2)
  )
  if(!is.null(id)) dt$Participant=id
  dt
}

rebc_to_circos_svs <- function(df, filter_statement=NA, sectors=NULL){
  
  if(is.null(sectors)) sectors <- c(paste0("chr", c(1:22, "X", "Y")))
  
  df <- df %>% ungroup()
  # thanks to https://stackoverflow.com/a/61692420/4783389
  if(!is.na(filter_statement)) df <- df %>% filter(eval(rlang::parse_expr(filter_statement)))
  df <- df %>%
    relocate(chr1, pos1, chr2, pos2, VAF) %>% # TODO update for hg38 chr
    mutate(chr1 = paste0("chr", chr1), chr2=paste0("chr", chr2)) %>%
    mutate(chr1 = case_when(chr1=="chr23" ~ "ChrX",
                            chr1=="chr24" ~ "ChrY",
                            TRUE ~ chr1)) %>%
    filter(chr1 %in% sectors & chr2 %in% sectors) %>%
    mutate(col=case_when( class == "deletion" ~ "#08519C", # blues
                          class == "tandem_dup" ~ "#006D2C", # green
                          class == "long_range" ~ "#54278F", # purple
                          class == "inversion" ~ "#006D2C", # green
                          class == "inter_chr" ~ "#3B3B3B")) # grey
  df
}

read_acs = function(acs, id=NULL){
  tsv <- read_tsv(acs)
  if(!is.null(id)) {
    tsv$Tumor_Sample_Barcode=id
    tsv <- tsv |> relocate(Tumor_Sample_Barcode)
  }
  tsv
}


merge_acs_pair = function(acs1, acs2, contigs=c(1:22, c("X", "Y")), id=NULL){
  acs <- bind_rows(acs1, acs2) |>
    filter(Chromosome %in% contigs) |>
    mutate(Chromosome = fct_recode(factor(Chromosome, levels=contigs))) |>
    arrange(Chromosome, Start.bp)
  if(!is.null(id)) {
    acs$Pair=id
    acs <- acs |> relocate(Pair)
  }
  acs
}


colorize.snv.per.classification <- function(type, ns=NULL) {
  if(is.null(ns)) ns <- c("COULD_NOT_DETERMINE", "DE_NOVO_START_IN_FRAME", "DE_NOVO_START_OUT_FRAME", 
                          "Frame_Shift_Del", "In_Frame_Del", "Missense_Mutation", "Nonsense_Mutation", 
                          "Nonstop_Mutation", "Splice_Site", "START_CODON_SNP")
  if(type %in% ns & length(type)==1) "#F16913FF"
  else "#4EB3D380" # 80 is for 50% transparency in HEXA
}


snv.to.circos <- function(x, sectors=NULL) {
  if(is.null(sectors)) sectors <- c(paste0("chr", c(1:22, "X", "Y")))
  x <- x |> 
    rowwise() |>
    transmute(chr=Chromosome, # or chr+Chromosome if refbuild hg19 
              start=Start_position, 
              end=End_position, 
              ccf_hat, 
              Mut=Hugo_Symbol,
              Class=Variant_Classification) |>
    filter(chr %in% sectors)
  if(nrow(x)>=1) {
    x <- x |> rowwise() |> mutate(Color=colorize.snv.per.classification(Class)) }
  else {
    x$Color=NA # trick for empty column to avoid colorize.snv.per.classification bugs
  }
  x |>
      arrange(Color, chr, start, end, ccf_hat)
}


tx.to.bed1 <- function(x) {
  x %>% 
    mutate(IG_CHR=paste0("chr", IG_CHR)) %>%
    dplyr::select(IG_CHR, IG_START, IG_END)
}
tx.to.bed2 <- function(x) {
  x %>% 
    mutate(PARTNER_CHR=paste0("chr", PARTNER_CHR)) %>%
    dplyr::select(PARTNER_CHR, PARTNER_START, PARTNER_END)
}


cnv.acs.to.circos <- function(df, sectors=NULL){
  # for alleliccapseg results
  if(is.null(sectors)) sectors <- c(paste0("chr", c(1:22, "X", "Y")))
  df <- df %>%
    mutate(col=col_fun_tau(Segment_Mean),
           chr=paste0("chr", Chromosome)) %>%
    filter(chr %in% sectors) %>%
    transmute(chr, start=Start, end=End, Segment_Mean, col)
  df
}

cnv.seg.to.circos <- function(df, sectors=NULL){
  # for gatk
  if(is.null(sectors)) sectors <- c(paste0("chr", c(1:22, "X", "Y")))
  df <- df %>%
    mutate(col_l2r=col_fun_l2r(mean_log2_ratio),
           chr=paste0("chr", chr)) %>%
    filter(chr %in% sectors) %>%
    transmute(chr, start, end, mean_log2_ratio, col_l2r)
  df
}

plotCancerCircos <- function(tx=NULL, cna=NULL, quick=TRUE, onco.tx=NULL, coding.snv=NULL, main=NULL, print.gene.names=TRUE){
  
  if(quick) { plotType=c('axis', 'labels') } else { plotType=c('axis', 'labels', 'ideogram') }
  
  circos.clear()
  par(lwd = 0.5)
  circos.par("cell.padding" = c(0, 0, 0, 0),
             "start.degree" = 90,
             "track.height" = 0.05,
             "canvas.xlim" = c(-1.3, 1.3),
             "canvas.ylim" = c(-1.3, 1.3))
  
  sectors.names <- c(1:22, "X", "Y") # as they appear in circos
  sectors <- c(paste0("chr", sectors.names)) # as in vcf output and ucsc reference
  
  circos.initializeWithIdeogram(species = "hg19", major.by = 1E9, axis.labels.cex = 1E-4, chromosome.index = sectors, plotType = plotType )
  # circos.initialize()
  # 
  if(print.gene.names==TRUE && !is.null(coding.snv) && is.data.frame(coding.snv) && nrow(coding.snv)>=1) {
    circos.genomicLabels(coding.snv, labels.column = 6, side = "outside", cex = 0.5, labels_height = 0.5, track.margin = c(0,0))
  }
  
  # circos.initializeWithIdeogram(species = "hg19", major.by = 1E9, axis.labels.cex = 1E-4, chromosome.index = sectors, plotType = c('axis', 'labels'))
  
  
  # circos.genomicIdeogram()
  
  if(!is.null(cna) && is.data.frame(cna) && nrow(cna)>=1) {
    circos.genomicTrackPlotRegion(cna, ylim = c(0, 5), panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col=value$col_l2r, border = NA, ...)
    }, bg.border = NA)
    circos.genomicTrackPlotRegion(cna, ylim = c(0, 5), panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col=value$col_call, border = NA, ...)
    }, bg.border = NA)
  }
  
  if(!is.null(coding.snv) && is.data.frame(coding.snv) && nrow(coding.snv)>=1) {
    circos.genomicTrack(coding.snv, numeric.column = 4, ylim = c(0,1),
                        panel.fun = function(region, value, ...) {
                          # numeric.column is automatically passed to `circos.genomicPoints()`
                          circos.genomicPoints(region, value, ...)
                        }) 
  }
  if(!is.null(tx) && is.data.frame(tx) && nrow(tx)>=1 ) {
    circos.genomicLink(tx[, c("chr1", "start1", "end1")] , tx[, c("chr2", "start2", "end2")], lwd = 1/2*10*tx$VAF)
  }
  if(!is.null(onco.tx) && is.data.frame(onco.tx) && nrow(onco.tx)>=1 ) {
    circos.genomicLink(onco.tx[, c("chr1", "start1", "end1")], onco.tx[, c("chr2", "start2", "end2")], col="darkred")
  }
  
  mtext(main)
}


plot.circos <- function() {
  circos.clear()
  
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  # par(mar = c(0, 0, 0, 0))
  circos.par("cell.padding" = c(0, 0, 0, 0), "start.degree" = 90, "track.height" = 0.05)
  # 
  # circos.initializeWithIdeogram(species = "hg19", 
  #                               major.by = 1E8, axis.labels.cex = 1E-4, 
  #                               chromosome.index = sectors, plotType = c('axis', 'labels'))
  # 
  circos.initializeWithIdeogram(species = "hg19", 
                                major.by = 5E7, axis.labels.cex = 1e-4, 
                                # chromosome.index = sectors, plotType = c('axis', 'ideogram', 'labels'))
                                chromosome.index = sectors, plotType = c('axis', 'labels'))
  
  # cn <- CTF.cnv.pairs.flat.annotated |> filter(participant == pair & tissue==itissue)
  # cn <- cnv.acs.to.circos(cn)
  
  CN.bin <- !is.null(cn) && is.data.frame(cn) && nrow(cn)>=1
  SV.rect <- !is.null(sv.rect) && is.data.frame(sv.rect) && nrow(sv.rect)>=1
  
  if(CN.bin & SV.rect) {
    bed_list = list(cn, sv.rect[,c("chr1", "pos1", "pos2", "VAF", "col")])
  } else if (CN.bin & !SV.rect) {
    bed_list = list(cn)
  } else if (!CN.bin & SV.rect){
    bed_list = list(sv.rect[,c("chr1", "pos1", "pos2", "VAF", "col")])
  } else {
    bed_list=NA
  }
  
  if(length(bed_list)>0) {
    circos.genomicTrackPlotRegion(bed_list[[1]], 
                                  # numeric.column = c(4, 4),
                                  ylim = c(0, 10), 
                                  panel.fun = function(region, value, ...) {
                                    # i = getI(...)
                                    circos.genomicRect(region, value, col=value$col, border = NA, ...)
                                    # circos.genomicRect(region, value, col=value$col, border = NA, ...)
                                  }, bg.border = NA)
  }
  
  # 
  #   if(!is.null(sv.rect) && is.data.frame(sv.rect) && nrow(sv.rect)>=1) {
  #     # circos.genomicTrackPlotRegion(cn, ylim = c(0, 10), panel.fun = function(region, value, ...) {circos.genomicRect(region, value, col=value$col_l2r, border = NA, ...)
  #     # }, bg.border = NA)
  #     circos.genomicTrackPlotRegion(sv.rect[,c("chr1", "pos1", "pos2", "col")], ylim = c(0, 10),
  #                                   panel.fun = function(region, value, ...) {
  #                                     circos.genomicRect(region, value, col=value$col, border = NA, ...)
  #                                   }, bg.border = NA)
  #   }
  #   
  
  if(!is.null(bulk.snv) && is.data.frame(bulk.snv) && nrow(bulk.snv)>=1) {
    circos.genomicTrack(bulk.snv, numeric.column = 4, ylim = c(0,1),
                        panel.fun = function(region, value, ...) {
                          # numeric.column is automatically passed to `circos.genomicPoints()`
                          circos.genomicPoints(region, value, cex = .3, pch = 16, col = value$Color, ...)
                        },bg.border = NA)
  }
  
  if(!is.null(snv) && is.data.frame(snv) && nrow(snv)>=1) {
    # circos.genomicTrack(snv, numeric.column = 4, ylim = c(0,1),
    #                     panel.fun = function(region, value, ...) {
    #                       # numeric.column is automatically passed to `circos.genomicPoints()`
    #                       circos.genomicPoints(region, value, ...)
    #                     },bg.border = NA)
    circos.genomicLabels(snv, labels.column = 5, side = "inside", cex = 0.5, labels_height = 0.5, track.margin = c(0,0))
  }
  
  
  if(!is.null(sv.link) && is.data.frame(sv.link) && nrow(sv.link)>=1 ) {
    circos.genomicLink(sv.link[, c("chr1", "pos1", "pos1")], 
                       sv.link[, c("chr2", "pos2", "pos2")],
                       col = sv.link$col,
                       rou = .7, lwd = .5)
  }
  
  
  if(!is.null(tx) && is.data.frame(tx) && nrow(tx)>=1 ) {
    circos.genomicLink(tx[, c("chr1", "start1", "end1")], tx[, c("chr2", "start2", "end2")], col =
                         "darkred", rou = .75, lwd = 2, h.ratio=.3)
  }
  
}

plot.ccf.density <- function(p.id, maf, 
                             maf.labels=TRUE, 
                             arrow=FALSE, 
                             save=TRUE, 
                             plot=TRUE, 
                             clustering="pre", 
                             genes.of.interest=NULL, 
                             x.suffix="BMPCs",
                             y.suffix="CMMCs",
                             x.label="Bone marrow",
                             y.label="Peripheral blood",
                             n.bins=30,
                             hex.color="darkblue",
                             min.COSMIC.total.alterations.in.gene=-1, 
                             non.silent.only=TRUE) {
  maf.pid <- subset(maf, Patient_ID==p.id)
  gplot.1 <- ggplot(maf.pid)
  if("pre" %in% clustering) {
    gplot.1 <- gplot.1 +
      geom_hex(bins=n.bins, 
               aes_string(x=paste0("preDP_ccf_mean_", x.suffix), y=paste0("preDP_ccf_mean_", y.suffix))) + 
      scale_fill_gradient(low="white", high = first(hex.color))
  }
  if ("post" %in% clustering) { # NOT WORKING
    warning("option post clustering does not work. Also, do not combine pre and post clustering currently.")
    gplot.1 <- gplot.1 +
      geom_hex(bins=n.bins, 
               aes_string(paste0("clust_ccf_mean_", x.suffix), paste0("clust_ccf_mean_", y.suffix))) + 
      scale_fill_gradient(low="white", high = last(hex.color))
  }
  if(maf.labels == TRUE) {
    maf.labels <- maf.pid
    if(!is.null(genes.of.interest)) maf.labels <- maf.labels %>% filter(Hugo_Symbol %in% genes.of.interest | COSMIC_total_alterations_in_gene > min.COSMIC.total.alterations.in.gene)
    if(non.silent.only) maf.labels <- maf.labels %>% filter(Variant_Classification %in% non.synonymous)
    message("Found ", nrow(maf.labels), " events to annotate out of ", nrow(maf.pid), " total variants.")
    if(nrow(maf.labels)>0) {
      maf.labels$label <- paste0(maf.labels$Hugo_Symbol, "\n", maf.labels$Protein_Change)
      if(arrow==TRUE) { #not working
        warning("option arrow TRUE is not working.")
        gplot.1 <- gplot.1 + 
          geom_text(data=maf.labels, 
                    aes(x+0.1, y+0.05, label=paste0(Hugo_Symbol, "\n", Protein_Change)),
                    size=2, hjust="left", vjust="middle") +
          geom_segment(data=maf.labels,
                       aes(xend=x, yend=y, x=x+0.1, y=y+0.05),
                       size=0.5, arrow = arrow(length = unit(0.03, "npc")))
      } else {
        # always post clustering for mutation annotation
        gplot.1 <- gplot.1 + 
          geom_point(data=maf.labels,
                     aes_string(x=paste0("clust_ccf_mean_", x.suffix), 
                                y=paste0("clust_ccf_mean_", y.suffix)), size=0.3) +
          geom_text(data = maf.labels, 
                    aes_string(paste0("clust_ccf_mean_", x.suffix), 
                               paste0("clust_ccf_mean_", y.suffix), 
                               label="label"),
                    position = position_jitter(width = .02, height = .02),
                    size=2, hjust="left", vjust="middle") 
      }
    }
  }  
  gplot.1 <- gplot.1 +  
    labs(x=x.label, y=y.label, fill="# SNVs") +
    scale_x_continuous(limits=c(0, NA), labels=scales::percent) +
    scale_y_continuous(limits=c(0, NA), labels=scales::percent) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          aspect.ratio = 1, 
          legend.position = "none")
  if(save==TRUE) {
    ggsave(paste0("figures/ccf-graphs_", p.id, "_scatter.png"), width = 3, height = 3)
    ggsave(paste0("figures/ccf-graphs_", p.id, "_scatter.pdf"), width = 3, height = 3)
  }
  if (plot==TRUE){
    gplot.1
  }
}



hrd.jb <- function(df){
  # cols 1 and 2 reserved for sample id and cohort
  # JB definition of HRD follows IFM trisomy of >= 2 chr
  # very sad to code like this but it's fast and it works
  df$amp_chr_3_HD <- ifelse(df$amp_chr_3p > 0.1 & df$amp_chr_3q > 0.1, 1, 0)
  df$amp_chr_5_HD <- ifelse(df$amp_chr_5p > 0.1 & df$amp_chr_5q > 0.1, 1, 0)
  df$amp_chr_7_HD <- ifelse(df$amp_chr_7p > 0.1 & df$amp_chr_7q > 0.1, 1, 0)
  df$amp_chr_9_HD <- ifelse(df$amp_chr_9p > 0.1 & df$amp_chr_9q > 0.1, 1, 0)
  df$amp_chr_11_HD <- ifelse(df$amp_chr_11p > 0.1 & df$amp_chr_11q > 0.1, 1, 0)
  df$amp_chr_15_HD <- ifelse(df$amp_chr_15q > 0.1, 1, 0)
  df$amp_chr_19_HD <- ifelse(df$amp_chr_19p > 0.1 & df$amp_chr_19q > 0.1, 1, 0)
  df$amp_chr_21_HD <- ifelse( df$amp_chr_21q > 0.1, 1, 0)
  df$total_HD_chr <- rowSums(df[,c("amp_chr_3_HD",  "amp_chr_5_HD", "amp_chr_7_HD", "amp_chr_9_HD", "amp_chr_11_HD", "amp_chr_15_HD", "amp_chr_19_HD", "amp_chr_21_HD")])
  df$HRD_JB <- df$total_HD_chr >= 2
  df$HRD_JB_values <- ifelse(df$HRD_JB, paste0("HRD(",
                                               ifelse(df$amp_chr_3_HD, "3,", ""), 
                                               ifelse(df$amp_chr_5_HD, "5,", ""), 
                                               ifelse(df$amp_chr_7_HD, "7,", ""), 
                                               ifelse(df$amp_chr_9_HD, "9,", ""), 
                                               ifelse(df$amp_chr_11_HD, "11,", ""), 
                                               ifelse(df$amp_chr_15_HD, "15,", ""), 
                                               ifelse(df$amp_chr_19_HD, "19,", ""), 
                                               ifelse(df$amp_chr_21_HD, "21,", ""), 
                                               ")"),"")
  df$HRD_JB_values <- str_replace(df$HRD_JB_values, pattern = ",\\)", replacement = "\\)") # so sorry for that horror
  df %>% select(1, 2, total_HD_chr, HRD_JB, HRD_JB_values)
}
