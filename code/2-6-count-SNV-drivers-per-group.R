setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

# NOTES FROM 2024-03-30 -------------------------
# RUN 2-first

final.snv.numeric.matrix[1:5, 1:5]
patient.annotation.df[1:5, 1:5]

# Distribution of MM driver point mutations from MGUS, SMM, to MM
# for now only in >1%

# count SNV drivers --------------

N.drivers.per.participant <- colSums(1*(final.snv.numeric.matrix>=1))

count.snv.drivers.per.participant <- data.frame(drivers=N.drivers.per.participant)  |> 
  rownames_to_column("ID") |>
  inner_join(patient.annotation.df|>rownames_to_column("ID")) |>
  group_by(drivers, Stage) |>
  count()

count.snv.drivers.per.participant

n.drivers.per.stage.plot <- ggplot(count.snv.drivers.per.participant, aes(drivers, n, fill=Stage)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks=0:20) +
  facet_grid(Stage~., scales = "free_y", switch = "both") +
  scale_fill_manual(values=brewer.pal(7, "Set1")) +
  labs(x="Number of MM drivers mutated", y="Frequency") +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave("../figures/n.drivers.per.stage.plot.pdf", width = 2, height = 2)
ggsave("../figures/n.drivers.per.stage.plot.png", width = 2, height = 2)