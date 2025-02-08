library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(RColorBrewer)

file <- "../../Alberge et al 2024 MMap_raw data archive/DocsSubmissionJuly2024/editor revisions jan 13 2025/AlbergeDuttaPoletti_SupplementaryTables.xlsx"
st1 <- read_xlsx(file, skip = 9)
head(st1)

table(st1$Cohort)

st1 <- st1 |>
  filter(Cohort %in% c("1-SU2C", "2-CTC")) |>
  mutate(Assay=case_when(Patient%in%c(paste0("SCI", str_pad(pad="0",side = "left",width = 3, 1:29)))~"Bone marrow\nCellSearch (R)",
                         Cohort=="1-SU2C"~"Bone marrow\nautoMACS (R)",
                         Cohort=="2-CTC"~"Peripheral Blood\nCellSearch (R)",
                         TRUE~"Unexpected"
                         ))

table(st1$Assay)
  
stats <- dunn_test(st1, `Purity (ABSOLUTE)`~Assay, p.adjust.method="bonferroni") %>% 
  add_xy_position(step.increase = 3) %>%
  mutate(p.adj=format(p.adj, digits=2))

summary.stats = st1 |> 
  group_by(Assay) |>
  summarize(n=n(), 
            median=median(`Purity (ABSOLUTE)`), 
            mad=mad(`Purity (ABSOLUTE)`), 
            text=paste0("N = ", n, "\n", "med. ", scales::percent(median), "\nmad ", scales::percent(mad)))

fig.supp.cs.facs <- ggplot(st1, aes(Assay)) +
  geom_boxplot(aes(y=`Purity (ABSOLUTE)`, color=Assay)) +
  geom_jitter(aes(y=`Purity (ABSOLUTE)`),shape=21,size=0.5, color="lightgrey", alpha=0.7) +
  geom_text(data=summary.stats, aes(label = text), y=0, size = 6, size.unit = "pt") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=brewer.pal(3, "Set2")) +
  scale_y_continuous(limits=c(0, NA), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent) +
  # stat_compare_means(y=aes(`Purity (ABSOLUTE)`)) +
  stat_pvalue_manual(stats, label = "p.adj", size = 2.7) +
  theme(text = element_text(size=6))

fig.supp.cs.facs

ggsave("figures/supp_fig_cs_facs_purity_comp.png", width = 2.5, height = 2.5)  
ggsave("figures/supp_fig_cs_facs_purity_comp.pdf", width = 2.5, height = 2.5)  
 