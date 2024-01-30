library(tidyverse)
library(pwr)

# Aim 1
d_CCEP_sig <- data.frame()
CCEP_data <- CCEP_amy_avg
timecourse_sig <- c(119:165,349:430)
for (i in 1:length(timecourse_sig)){
  print(i)
  CCEP_timepoint <- filter(CCEP_data, Time == timecourse_sig[i])
  mean <- mean(CCEP_timepoint$CCEP_norm)
  sd <- sd(CCEP_timepoint$CCEP_norm)
  mu <- 0
  es <- data.frame(abs(mean-mu)/sd)
  d_CCEP_sig <- bind_rows(d_CCEP_sig,es)
}
mean(d_CCEP_sig$abs.mean...mu..sd)
range(d_CCEP_sig$abs.mean...mu..sd)

# Aim 2a
TEP_data <- TEP_amy_avg
timecourse_sig <- c(52:91,194:244,326:347)
d_TEP_sig <- data.frame()
for (i in 1:length(timecourse_sig)){
  print(i)
  TEP_tp_tms <- filter(TEP_data, Time == timecourse_sig[i] & Condition == "tms")
  TEP_tp_sham <- filter(TEP_data, Time == timecourse_sig[i] & Condition == "sham")
  mean_tms <- mean(TEP_tp_tms$iTEP_norm)
  mean_sham <- mean(TEP_tp_sham$iTEP_norm)
  sd_tms <- sd(TEP_tp_tms$iTEP_norm)
  sd_sham <- sd(TEP_tp_sham$iTEP_norm)
  sd <- sqrt((sd_tms^2+sd_sham^2)/2)
  es <- data.frame(abs(mean_tms-mean_sham)/sd)
  d_TEP_sig <- bind_rows(d_TEP_sig,es)
}
mean(d_TEP_sig$abs.mean_tms...mean_sham..sd)
range(d_TEP_sig$abs.mean_tms...mean_sham..sd)


# Aim 2b
TMSres_rM1_AMY <- filter(TMSres[c(3,4,5,7,8)], site == "R_M1")
mean_lpDLPFC <- mean(TMSres_lpDLPFC_AMY$FIRST_B_amyg_small)
mean_rM1 <- mean(TMSres_rM1_AMY$FIRST_B_amyg_small)
sd_lpDLPFC <- sd(TMSres_lpDLPFC_AMY$FIRST_B_amyg_small)
sd_rM1 <- sd(TMSres_rM1_AMY$FIRST_B_amyg_small)
sd <- sqrt((sd_lpDLPFC^2+sd_rM1^2)/2)
es <- abs(mean_lpDLPFC-mean_rM1)/sd
es

# Aim 3
pwr.r.test(r= 0.54, sig.level = 0.05, power = 0.8)
pwr.r.test(r= 0.25, sig.level = 0.05, power = 0.8)

pwr.r.test(r= 0.54, sig.level = 0.05, n = 21)
pwr.r.test(r= 0.25, sig.level = 0.05, n = 69)

pwr.r.test(r= 0.44, sig.level = 0.05, power = 0.8)
pwr.r.test(r= 0.28, sig.level = 0.05, power = 0.8)
