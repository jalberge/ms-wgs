library(tidyverse)
library(RColorBrewer)

TOTAL.DRAWS = 100
TOTAL.N = 150

PRE.N = c(30, 60)

PRE.CCF = seq(0, 1, by=0.01)
POST.CCF = seq(0, 1, by=0.01)

DATA = expand_grid(DRAW=1:TOTAL.DRAWS, TOTAL.N=TOTAL.N, PRE.N=PRE.N, POST.N=POST.N, PRE.CCF=PRE.CCF, POST.CCF=POST.CCF)

RANDOM.DATA <- DATA |>
  rowwise() |>
  mutate(POST.N = TOTAL.N-PRE.N,
         GROUP = paste0("Case: ", PRE.N, "/",  POST.N, " Cells"),
         PRE.DRAW=rbinom(1, PRE.N, PRE.CCF),
         POST.DRAW=rbinom(1, POST.N, POST.CCF))

POWER <- RANDOM.DATA |> 
  mutate(
    P.val.less = fisher.test(matrix(c(PRE.DRAW, PRE.N-PRE.DRAW, POST.DRAW, POST.N-POST.DRAW), nrow = 2, byrow = TRUE), alternative = "l")$p.value,
    P.val.more = fisher.test(matrix(c(PRE.DRAW, PRE.N-PRE.DRAW, POST.DRAW, POST.N-POST.DRAW), nrow = 2, byrow = TRUE), alternative = "g")$p.value) |>
  group_by(TOTAL.N, PRE.N, POST.N, PRE.CCF, POST.CCF, GROUP) |>
  summarise(Power.less = sum(P.val.less<0.05)/TOTAL.DRAWS,
            Power.more = sum(P.val.more<0.05)/TOTAL.DRAWS)

LINES <- POWER |> 
  group_by(PRE.CCF, GROUP) |>
  summarise(
    MIN.POST.CCF.80.POWER.LESS = min(1, POST.CCF[Power.less > 0.8]),
    MIN.POST.CCF.80.POWER.MORE = max(0, POST.CCF[Power.more > 0.8]))

LINES <- LINES |> 
  mutate(GROUP=case_when(GROUP=="Case: 30/120 Cells" ~ "Case 2: 30/120 Cells",
                         GROUP=="Case: 60/90 Cells" ~ "Case 1: 60/90 Cells"))

# ggplot(LINES |> filter(GROUP=="Pre: 30; Post: 120 Cells"), aes(PRE.CCF,  fill=GROUP)) +
ggplot(LINES, aes(PRE.CCF, fill=GROUP)) +
  geom_ribbon(aes(ymin=MIN.POST.CCF.80.POWER.LESS, ymax=1), position = "identity") +
  geom_ribbon(aes(ymin=0, ymax=MIN.POST.CCF.80.POWER.MORE), position = "identity") +
  geom_segment(aes(x=0.05, y=0.05, xend=0.06, yend=0.26)) +
  scale_x_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_y_continuous(limits=c(0, 1), labels = scales::percent) + 
  scale_fill_manual(values=brewer.pal(8, "Blues")[c(4, 7)]) +
  labs(x="Fraction of cells\nbefore progression", 
       y="Fraction of cells\nafter progression",
       fill="Area powered (80%) for detection") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "top")

ggplot(POWER, aes(PRE.CCF, POST.CCF, fill=Power)) +
  geom_raster() +
  # geom_point() +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio = 1)


