# run 2-2-redo-score


# Load surv ---------------------------------------------------------------


comm.surv <- read_tsv("../../../MMRF_IA20/CoMMpass_IA20_FlatFiles/MMRF_CoMMpass_IA20_STAND_ALONE_SURVIVAL.tsv")

mini.comm.surv <- comm.surv |>
  select(PUBLIC_ID, any_of(c( "deathdy", "lstalive", "lvisitdy", "linesdy1",
                              "pfsdy", "ttpfs", "censpfs", "pfscdy", "ttcpfs", "ttpfsw", 
                             "ttcpfsw", "pfsdy1", "ttpfs1", "censpfs1", "pfs1cdy", "ttcpfs1", 
                             "ttpfs1w", "ttcpfs1w", "ttcos", "censos")))

head(mini.comm.surv)

# first line: ttcpfs1, 1-censpfs1


# Merge with score --------------------------------------------------------


# add mm-like score
sig.comm <- sig.mm |> filter(startsWith(ID, "MMRF"))

mini.comm.surv.sig <- inner_join( mini.comm.surv, sig.comm, by=c("PUBLIC_ID"="ID"))
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.33, 0.67, 1)), 
                                    include.lowest = TRUE)
table(mini.comm.surv.sig$sig_class, mini.comm.surv.sig$total_sig, useNA = 'ifany')
# coxph(Surv(ttcpfs1, censpfs1) ~ total_sig, data=mini.comm.surv.sig)


# Scores PFS ------------------------------------------------------------------


## quartile --------------------------------------------------------------


# score quartile
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)

## tertile --------------------------------------------------------------

# score tertile
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.33, 0.67, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)
# HR 1.76 pour score >= 5, 1.16 pour >=3 <=4

## median --------------------------------------------------------------

# score median
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, 
                                    breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.5, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
summary(coxph(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig))
# HR 1.5

# comp with MAF for instance
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ t14_16, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcpfs1, censpfs1) ~ t14_16, data=mini.comm.surv.sig)
# absolutely no signal HR ~1

# comp with 1q
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ amp_chr_1q, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcpfs1, censpfs1) ~ amp_chr_1q, data=mini.comm.surv.sig)
# good signal HR 1.32

# comp with MYC
fit <- survfit(Surv(ttcpfs1, censpfs1) ~ amp_chr_8q, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcpfs1, censpfs1) ~ amp_chr_8q, data=mini.comm.surv.sig)
# good signal HR 1.31




# Scores OS ------------------------------------------------------------------


## quartile --------------------------------------------------------------


# score quartile
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
# HR 2.0459 pour le haut groupe >= 5

## tertile --------------------------------------------------------------

# score tertile
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.33, 0.67, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
# HR 2.0460 pour score >= 5

## median --------------------------------------------------------------

# score median
mini.comm.surv.sig$sig_class <- cut(mini.comm.surv.sig$total_sig, 
                                    breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.5, 1)), 
                                    include.lowest = TRUE)
fit <- survfit(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
summary(coxph(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig))
# HR 1.6555 pour >=4

# comp with MAF for instance
fit <- survfit(Surv(ttcos, censos) ~ t14_16, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcos, censos) ~ t14_16, data=mini.comm.surv.sig)
# absolutely no signal HR ~1

# comp with 1q
fit <- survfit(Surv(ttcos, censos) ~ amp_chr_1q, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcos, censos) ~ amp_chr_1q, data=mini.comm.surv.sig)
# good signal HR 1.6

# comp with MYC
fit <- survfit(Surv(ttcos, censos) ~ amp_chr_8q, data=mini.comm.surv.sig)
ggsurvplot(fit, data=mini.comm.surv.sig)
coxph(Surv(ttcos, censos) ~ amp_chr_8q, data=mini.comm.surv.sig)
# good signal HR 1.2


# OS / PFS figure tertile -------------------------------------------------
legend.labs=c("Low", "Int", "High")
mini.comm.surv.sig$sig_class <- factor(cut(mini.comm.surv.sig$total_sig, breaks = quantile(mini.comm.surv.sig$total_sig, probs = c(0, 0.33, 0.67, 1)), 
                                    include.lowest = TRUE, labels = FALSE), labels=legend.labs)

## PFS------
fit1 <- survfit(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)
pfs.plot <- ggsurvplot(fit1, data=mini.comm.surv.sig, 
                       xlab="Time to PFS (years)",
                       ylab="Progression-free survival (PFS)",
                       xscale="d_y", 
                       break.time.by=365.25,
                       pval = TRUE,
                       surv.median.line="hv",
                       palette="grey",
                       legend.labs=legend.labs,
                       legend.title="MM-like score")
pfs.plot
coxph(Surv(ttcpfs1, censpfs1) ~ sig_class, data=mini.comm.surv.sig)


fit2 <- survfit(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
osplot <- ggsurvplot(fit2, data=mini.comm.surv.sig, 
                     xlab="Time to OS (years)",
                     ylab="Overall survival (PFS)",
                     xscale="d_y", 
                     break.time.by=365.25,
                     pval = TRUE,
                     surv.median.line="hv",
                     palette="grey",
                     legend.labs=legend.labs,
                     legend.title="MM-like score")
osplot

coxph(Surv(ttcos, censos) ~ sig_class, data=mini.comm.surv.sig)
# HR 2.0460 pour score >= 


library(patchwork)

os.pfs.plots <- pfs.plot$plot + osplot$plot

pdf("../figures/pfs-os-commpass-mm-like-score.pdf", width = 8, height = 5)
os.pfs.plots
dev.off()

ggsave("../figures/pfs-os-commpass-mm-like-score.png", os.pfs.plots,dpi = 320, width = 8, height = 5)
