setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(rstatix)
library(RColorBrewer)

smm <- read_tsv("../data/export_MM-like_score_SMM.tsv")

smm$type <- ifelse(smm$cohort=="Bustoros", "WXS", "WGS")

table(smm$type)

# higher risk patients in the Bustoros paper (which we knew)
ggplot(smm, aes(type, total_sig)) +
  geom_boxplot()

# what is it driven by?
count.tables <- smm |>
  pivot_longer(cols = matches("^(t1|KRAS|NRAS|FAM|Hyper|del_|amp_)"), names_to = "mutation", values_to = "status") |>
  group_by(mutation, status, type) |>
  count()
fishers <- count.tables |>
  group_by(mutation) |>
  
  do(tidy(fisher.test( xtabs(n ~ status + type, data = .data) ), data = .x)) |>
  
  ungroup() |>
  mutate(q.value=paste0("q=", format(p.adjust(p.value, method="fdr"), digits=1))) |>
  
  mutate(mutation=str_replace(mutation, "amp_chr_", "+")) |>
  mutate(mutation=str_replace(mutation, "del_chr_", "-")) |>
  mutate(mutation=str_replace(mutation, "t14_16", "t(14;16)(MAF)")) |>
  mutate(mutation=str_replace(mutation, "t14_20", "t(14;20)(MAFB)")) 

confints <- count.tables |>
  mutate(status=abs(status)) |> # -1 for MAF
  pivot_wider(id_cols = c(mutation, type), names_from =  c(status), values_from = n) |>
  mutate_all(~{replace_na(., 0)}) |>
  
  mutate(tidy(binom.test(`1`, `1`+`0`))) |>
  
  mutate(mutation=str_replace(mutation, "amp_chr_", "+")) |>
  mutate(mutation=str_replace(mutation, "del_chr_", "-")) |>
mutate(mutation=str_replace(mutation, "t14_16", "t(14;16)(MAF)")) |>
  mutate(mutation=str_replace(mutation, "t14_20", "t(14;20)(MAFB)")) 
head(confints)

wgs.wxs.comparison <- ggplot(confints, aes(type, `1`/(`1`+`0`))) +
  geom_bar(stat="identity", aes(fill=type)) +
  geom_linerange(aes(ymin=conf.low, ymax=conf.high)) +
  geom_text(data=fishers, aes(x=1.5, y=0.9, label = q.value), size = 2.5) +
  facet_wrap(~mutation) +
  scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
  scale_fill_manual(values=brewer.pal(8, "PRGn")[c(2, 7)]) +
  labs(x="", y="Frequency of Event (95% CI)") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave("../figures/comparison.frequency.events.MMlike.WGS.WXS.png", wgs.wxs.comparison, width = 4, height = 3)
ggsave("../figures/comparison.frequency.events.MMlike.WGS.WXS.pdf", wgs.wxs.comparison, width = 4, height = 3)


