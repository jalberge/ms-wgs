library(ComplexHeatmap)
library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/Projects/ms-wgs/")

# mtx <- read_table("data/HG38VCFs_Nov8.matrix") %>% column_to_rownames("FILE") %>% as.matrix()
mtx <- read_table("data/ALL20230110.matrix") %>% column_to_rownames("FILE") %>% as.matrix()
clean.name <- function(x) return(str_trim(str_remove(basename(x), ".fingerprinted.vcf")))
rownames(mtx) <- sapply(rownames(mtx), clean.name)
colnames(mtx) <- sapply(colnames(mtx), clean.name)
mtx

rownames(mtx) <- sapply(rownames(mtx), function(x) { strsplit(x, split = "::")[[1]]})[1,]
colnames(mtx) <- sapply(colnames(mtx), function(x) { strsplit(x, split = "::")[[1]]})[1,]

pdf("figures/matrix.fingerprint.pdf", paper = "USr")
Heatmap(mtx, show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 3))
dev.off()

tibble(mtx)
