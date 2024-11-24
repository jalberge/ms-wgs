library(tidyverse)
library(RColorBrewer)
library(scales)
# from lawrence 2014
# binom.test( rbinom(1, N, p0), N, p = 5E-6, alternative = "l")

calculate_p0 <- function(mu=1.4E-6, f_g=3.9, L=1500) {
  # mu 1.3E-6 for myeloma # background mutation frequency in mutations per base
  # f_g = 3.9 for the 90th percentile # gene-specific mutation rate factor calculated by MutSigCV
  # L = 1500 gene length in coding base
  p0 = 1 - (1 - mu * f_g) ^ (3*L/4)
  p0
}

calculate_p1 <- function(p0, r=0.1, m=0.1) {
  # r frequency of non silent mutations in the population above background
  # m misdetection rate
  p1 = p0 + r*(1-m)
  p1
}

find_min_patients_p0 <- function(N, p0, P=5E-6){
  # k number of trials (patients) with a non silent mutation, 
  # P genome wide significance
  k <- 0
  p_val <- 1
  while(p_val > P) {
    if(k+1>N)(return(N))
    p_val <-  binom.test(k+1, N, p0, alternative = "g")$p.value
    k <- k + 1
  }
  return(k)
}

calculate_power_p1 <- function(min.patients, N, p1){
  pbinom(min.patients-1, N, p1, lower.tail = FALSE)
}

N=300
p0 <- calculate_p0()
p1 <- calculate_p1(p0, r = 0.05)
min.patients <- find_min_patients_p0(N, p0)
power <- calculate_power_p1(min.patients, N, p1)

# curves

power.to.detect.drivers <- expand_grid(N = seq(50, 5000, by=50), 
                                       r = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1)) %>%
  mutate(p0=calculate_p0(mu=1.4E-6, f_g=3.9, L=1500)) %>%
  rowwise() %>%
  mutate(
    p1 = calculate_p1(p0, r = r),
    min.patients = find_min_patients_p0(N, p0),
    power = calculate_power_p1(min.patients, N, p1)
  )


ggplot(power.to.detect.drivers, aes(N, power, color=factor(r, labels=percent(unique(r))))) + 
  geom_line() +
  geom_hline(yintercept = 0.9, color="darkgrey", linetype=3) +
  # geom_vline(xintercept = c(224, 1034), color="darkgrey", linetype=3) +
  scale_color_manual(values=brewer.pal(9, "Blues")[3:9]) +
  scale_x_log10() +
  labs(x="Number of TN pairs", y="Power for 90% candidate drivers", color="Fraction") +
  theme_bw()

binom.test(min.patients, N, p1, alternative = "l")



binom.test( rbinom(1, N, p1), N, p = 5E-6, alternative = "l")

binom.test( sum ( rbinom(N, 1, p1) ), N, p = 5E-6, alternative = "g")

binom.test()

rbinom(200, 1, p1)
