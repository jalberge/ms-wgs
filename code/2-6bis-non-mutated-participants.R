setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)


# What about non-mutated people -------------------------------------------

non.mutated.patients <- colnames(final.numeric.matrix[ , colSums(final.numeric.matrix[snvs,])==0 ])

do.not.order.columns.based.on.2 <- c("Trisomies", "t(MYC;IgH/K/L)", snvs)
all.drivers.to.order.on.2 <- rownames(gene.q.val[rownames(final.numeric.matrix)[-which(rownames(final.numeric.matrix) %in% do.not.order.columns.based.on)],])

variable.object.to.order.on.2 <- c("Stage", all.drivers.to.order.on.2)
numeric.matrix.to.order.samples.2 <- (-1*(final.numeric.matrix)) |> t() |> as.data.frame() |> rownames_to_column("Sample ID")

new.order.2 <- meta.df[non.mutated.patients,] |> 
  rownames_to_column("Sample ID") |>
  left_join(numeric.matrix.to.order.samples) |>
  arrange(!!! rlang::syms(variable.object.to.order.on)) |>
  pull(`Sample ID`)
# > new.order[new.order %nin% colnames(final.class.matrix)]


driver.sum.ha.non.muated = HeatmapAnnotation(foo = anno_barplot(matrix(nc=2, cbind(mm.drivers, driver.sum-mm.drivers)[new.order.2,]), 
                                                                gp = gpar(fill = c("#6A51A3",  "#BCBDDC"),
                                                                          col = NA),
                                                                bar_width = 1,
                                                                axis_param=list(gp=gpar(fontsize = 6))),
                                             
                                             height = unit(0.8, "cm"), 
                                             annotation_label = "Number of\ndrivers", 
                                             annotation_name_side = "left", 
                                             annotation_name_gp = gpar(fontsize = 6)
)


Hm.non.mutant <- Heatmap(final.class.matrix[which(rownames(final.numeric.matrix) %nin% snvs), new.order.2], 
                         # column_order = unlist(hm_col_order, use.names = FALSE),
                         column_order = new.order.2,
                         show_heatmap_legend = FALSE,
                         column_split=meta.df[new.order.2, c("IMWG")],
                         column_title_gp = gpar(fontsize = 7),
                         row_split=gene.q.val[rownames(final.numeric.matrix)[which(rownames(final.numeric.matrix) %nin% snvs)], "group"],
                         row_title_gp = gpar(fontsize = 7),
                         show_column_names = FALSE,
                         cluster_rows = FALSE,
                         cluster_column_slices = FALSE,
                         show_column_dend = FALSE,
                         # col = c("white", "black"),
                         col = class.values.colors,
                         # row_names_gp = gpar(fontsize = 5, col=row.names.colors, fontface=row.names.fontface),
                         row_names_gp = gpar(fontsize = 6, col="black", fontface=row.names.fontface),
                         border = TRUE,
                         # left_annotation = gene.q.ha,
                         # right_annotation = or.col, 
                         # bottom_annotation = paste0("N=",length(new.order.2))
                         top_annotation = driver.sum.ha.non.muated
)

# final.numeric.matrix[1:5, column_order(Hm)$Unclassified] %>% clipr::write_clip()
patient.annotation.non.mutated <- HeatmapAnnotation(df=patient.annotation.df[new.order.2, c("IMWG", "Stage")], 
                                                    col = list(IMWG=imwg.colors, Stage=stage.colors, Cohort=cohort.colors, Age=age.colors, Assay=assay.colors, Gender=gender.colors),
                                                    height = unit(1.5, "cm"),
                                                    simple_anno_size_adjust = TRUE,
                                                    annotation_name_gp= gpar(fontsize = 9))

ht.list.non.mutant <- patient.annotation.non.mutated %v% Hm.non.mutant

draw(ht.list.non.mutant)

patient.annotation.df[new.order.2, c("IMWG", "Stage")] |> group_by(Stage) |> count()

pdf("../figures/FIG_S1_NoMutations_Heatmap_v3.pdf", width = 5, height = 6)
draw(ht.list.non.mutant)
dev.off()
