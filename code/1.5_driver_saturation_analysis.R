setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# fit logistic regression models to observed counts of drivers (per group of frequency) 
# to test for plateaux with around 1,030 and 2% prevalence
library(tidyverse)
library(RColorBrewer)
library(xgxr)

# UTILS -------------------------------------------------------------------

set.seed(1)

predict.F <- function(model_fit, data){
  stats:::predict(model_fit, se.fit = TRUE, interval="confidence")
}

predict.T <- function(model_fit, data){
  xgxr:::predict.nls(model_fit, se.fit = TRUE, interval="confidence")
}


# SIGNATURE ANALYZER ------------------------------------------------------

l=list.files(path="../data/MUTSIG2CV/", pattern = "*summary_saturation.txt", full.names = TRUE)
saturation.analysis.multiple.runs <- do.call(rbind, lapply(l, read_tsv))
saturation.analysis <- saturation.analysis.multiple.runs %>%
  group_by(n_patients) %>%
  slice_sample(n=1)

long.saturation.analysis <- saturation.analysis|>
  pivot_longer(cols=-c(index, n_patients, total_drivers, list_drivers), 
               names_to = "Frequency", 
               values_to = "Drivers") |>
  mutate(Frequency=fct_recode(Frequency, 
                              "<1%"="n_below_1", 
                              "1-2%"="n_over_1", 
                              "2-3%"="n_over_2", 
                              "3-4%"="n_over_3", 
                              "4-5%"="n_over_4",
                              ">=5%"="n_over_5")) |>
  mutate(Drivers.noised=Drivers+rnorm(length(Drivers), mean = 0, sd = 0.3))

# tips
# data is new column from nesting
# predict nls doesn't have confidence intevrals
# see https://stackoverflow.com/questions/61341287/how-to-calculate-confidence-intervals-for-nonlinear-least-squares-in-r
# https://stackoverflow.com/questions/60663962/r-making-predictions-and-confidence-intervals-with-different-models-for-each-g
# https://stackoverflow.com/questions/65431315/multiple-linear-model-prediction-in-dplyr

nls.models.predictions <- long.saturation.analysis |>
  group_by(Frequency) |>
  nest() |> 
  mutate(model = data |> map( ~nls(Drivers.noised ~ SSlogis(n_patients, Asym, xmid, scal), data=.) )) |>
  # mutate(pred = map2(model, data, ~as.data.frame(predFit(.x, .y, interval = "confidence", type = "response", level=0.95)))) |>
  mutate(pred = map(model, list(predict.T, as.data.frame))) |>
  unnest(c(pred, data))

head(nls.models.predictions)

ribbon.palette <- brewer.pal(n=8, name = "Set1")[c(7, 5:1)]

mutsig.plot <- ggplot(data = nls.models.predictions, mapping = aes(n_patients, color=Frequency)) + 
  geom_jitter(aes(y=Drivers), width = 0, height = 0.5, alpha=0.2, size=0.5) +
  geom_ribbon(aes(ymin=fit.lwr, ymax=fit.upr, fill=Frequency), alpha=0.1, colour=NA) +
  geom_line(aes(y=fit.fit))+
  geom_line(aes(y=fit.lwr), alpha=0.3)+
  geom_line(aes(y=fit.upr), alpha=0.3)+
  labs(x="Cohort size (n)", y="Significant drivers\n(MutSig, q<0.1)") +
  scale_color_manual(values=ribbon.palette) +
  scale_fill_manual(values=ribbon.palette) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.7), panel.grid = element_blank())

mutsig.plot
ggsave(plot=mutsig.plot, "../figures/figures_mutsig_saturation.png", width = 4, height = 3, scale = 1.2)

# GISTIC ------------------------------------------------------------------

gistic.basic <- read_table("../data/agg_basic_stat_gistic.tsv")

# first when no driver at all n_patients is not reported
gistic.basic <- gistic.basic |> mutate(n_patients = case_when(is.na(n_patients)~index+11, TRUE~n_patients))

# second n_over_1 is na
gistic.basic <- gistic.basic |> mutate(across(starts_with("n_"), ~ pmax(0, .x, na.rm = TRUE) ))

set.seed(1)
long.gistic.analysis <- gistic.basic|>
  pivot_longer(cols=-c(index, n_patients, total_drivers, list_drivers), 
               names_to = "Frequency", 
               values_to = "Drivers") |>
  mutate(Frequency=fct_recode(Frequency, 
                              "<1%"="n_below_1", 
                              "1-2%"="n_over_1", 
                              "2-3%"="n_over_2", 
                              "3-4%"="n_over_3", 
                              "4-5%"="n_over_4",
                              ">=5%"="n_over_5")) |>
  mutate(Drivers.noised=Drivers+rnorm(length(Drivers), mean = 0, sd = .3))

# GISTIC doesn't have rare events therefore doesn't converge :(
# control=list(warnOnly=TRUE)

nls.models.gistic.predictions <- long.gistic.analysis |>
  group_by(Frequency) |>
  nest() |> 
  mutate(model = data |> map( ~nls(Drivers.noised ~ SSlogis(n_patients, Asym, xmid, scal), data=.))) |>
  # mutate(pred = map2(model, data, ~as.data.frame(predFit(.x, .y, interval = "confidence", type = "response", level=0.95)))) |>
  mutate(pred = map(model, list(predict.T, as.data.frame))) |>
  unnest(c(pred, data))

head(nls.models.gistic.predictions)

ribbon.palette <- brewer.pal(n=8, name = "Set1")[c(7, 5:1)]

gistic.plot <- ggplot(data = nls.models.gistic.predictions, mapping = aes(n_patients, color=Frequency)) + 
  geom_jitter(aes(y=Drivers), width = 0, height = 0.5, alpha=0.2, size=0.5) +
  geom_ribbon(aes(ymin=fit.lwr, ymax=fit.upr, fill=Frequency), alpha=0.1, colour=NA) +
  geom_line(aes(y=fit.fit))+
  geom_line(aes(y=fit.lwr), alpha=0.3)+
  geom_line(aes(y=fit.upr), alpha=0.3)+
  labs(x="Cohort size (n)", y="Significant drivers\n(GISTIC2, q<0.1)") +
  scale_color_manual(values=ribbon.palette) +
  scale_fill_manual(values=ribbon.palette) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.7), panel.grid = element_blank())

gistic.plot

ggsave(plot=gistic.plot, "../figures/figures_gistic_saturation.png", width = 4, height = 3, scale = 1.2)
