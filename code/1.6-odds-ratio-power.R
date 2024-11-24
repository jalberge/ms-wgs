setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(epiR)
library(tidyverse)
library(RColorBrewer)

total.n.drivers <- 100 # rounding from 103

# see https://cran.r-project.org/web/packages/epiR/vignettes/epiR_sample_size.html
epi.sscohortc(N = NA, # eligible (infinite)
              irexp1 = 2*0.15, # incidence risk in exposed
              irexp0 = 0.15,# incidence risk in unexposed
              power = NA, #power
              n = 1030, # total number of individuals
              r = 812/218, # exposed / unexposed
              design = 1, 
              sided.test = 2, 
              nfractional = FALSE, 
              conf.level = 1-.05/total.n.drivers)$power

input.space <- expand_grid(irexp0 = seq(0.01, 0.5, by=0.01),
                           irr = c(0.5, 1.5, 2, 2.5))

power.results <- input.space |>
  mutate(index=row_number()) |>
  group_by(index, irexp0) |>
  group_modify(~{data.frame(epi.sscohortc(N = NA, # eligible (infinite)
                                irexp1 = .x$irr*.x$irexp0, # incidence risk in exposed
                                irexp0 = .x$irexp0,# incidence risk in unexposed
                                power = NA, #power
                                n = 1030, # total number of individuals
                                r = 812/218, # exposed / unexposed
                                design = 1, 
                                sided.test = 2, 
                                nfractional = FALSE, 
                                conf.level = 1-.05/total.n.drivers))}, .keep = TRUE)

fig1 <- ggplot(power.results, aes(irexp0, power, color=factor(irr))) +
  
  geom_hline(yintercept = 0.8, linetype=2, color="darkgrey") +
  
  geom_point(size=0.5, shape=3) +
  geom_line(size=0.5) +
  
  scale_x_log10(labels=scales::percent) +
  scale_color_manual(values=rev(brewer.pal(11, "RdBu")[c(2, 4, 8, 10)])) +
  
  labs(x="Incidence of mutation in MGUS/SMM", 
       y=str_c("Detection power at α = 0.05/100"),
       color="Incidence risk ratio") +
  
  theme_bw(base_size = 7) +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.2, 0.7), 
        legend.title = element_text(size=6), 
        legend.text = element_text(size=6),
        legend.background = element_rect(size=0.25, color="#333333"))

fig1 

ggsave("../figures/revisions_power_to_detect_mutated_genes_odds_ratio.png", width = 3, height = 3)  
