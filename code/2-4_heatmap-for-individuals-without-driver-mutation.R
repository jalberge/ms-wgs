setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

# NOTES FROM 2024-03-30 -------------------------
# RUN 2.1 first

# What about non-mutated individuals --------------------------------------

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

Hm.non.mutant <- Heatmap(final.class.matrix[which(rownames(final.numeric.matrix) %nin% snvs),new.order.2], 
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
                         col = class.values.colors,
                         row_names_gp = gpar(fontsize = 5, col="black", fontface=row.names.fontface),
                         border = TRUE
)

patient.annotation.non.mutated <- HeatmapAnnotation(df=patient.annotation.df[new.order.2, c("IMWG", "Stage")], 
                                                    col = list(IMWG=imwg.colors, Stage=stage.colors, Cohort=cohort.colors, Age=age.colors, Assay=assay.colors, Gender=gender.colors),
                                                    height = unit(1.5, "cm"),
                                                    simple_anno_size_adjust = TRUE,
                                                    annotation_name_gp= gpar(fontsize = 9))

ht.list.non.mutant <- patient.annotation.non.mutated %v% Hm.non.mutant

draw(ht.list.non.mutant)

patient.annotation.df[new.order.2, c("IMWG", "Stage")] |> group_by(Stage) |> count()

pdf("../figures/FIG_S1_NoMutations_Heatmap_v2.pdf", width = 5, height = 5)
draw(ht.list.non.mutant)
dev.off()
