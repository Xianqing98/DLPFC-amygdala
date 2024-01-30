library(tidyverse)
library(readxl)

# Prepare dataset that suit EEGLAB input
# It should be in a [Channel x Timepoint] matrix, un-epoched

rootpath <- "/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/"
wd <- paste0(rootpath,"RData_CCEP/CCEP_csv_L-DLPFC_amy") ## unfiltered data
setwd(wd)

## Read patients csv
CCEP_amy_LDLPFC_orig <- read_csv(list.files(wd, pattern = '*.csv'))
CCEP_amy_LDLPFC_orig <- filter(CCEP_amy_LDLPFC_orig, Time <= 1500 & Time > -1000)
CCEP_amy_LDLPFC_orig$Channel <- paste(CCEP_amy_LDLPFC_orig$Patient, CCEP_amy_LDLPFC_orig$Channel, sep = "_")

## Add channel info
chlinfo_amy <- read_excel(paste0(rootpath, "Selected_Channel_AMY.xlsx"),
                          sheet = "esTT L-DLPFC") # l dlpfc
chlinfo_amy$Channel <- paste(chlinfo_amy$Patient, chlinfo_amy$Channel, sep = "_LFPx")
CCEP_amy_LDLPFC_orig <- left_join(CCEP_amy_LDLPFC_orig, chlinfo_amy[c(2,14,22)])

## delete questionable channels
### L-DLPFC
CCEP_amy_LDLPFC_orig <- filter(CCEP_amy_LDLPFC_orig, 
                               Channel != "634_LFPx55" & 
                               Channel != "634_LFPx56" )
ChannelList <- CCEP_amy_LDLPFC_orig %>% 
  group_by(Channel) %>% summarise(Patient = mean(Patient))
ChannelList <- ChannelList[[1]]

## Rearrange the iEEG data to wide format
CCEP_amy_LDLPFC <- pivot_wider(CCEP_amy_LDLPFC_orig, 
                               names_from = c(Time, TrialNumber), 
                               values_from = Amplitude)

## Save as TF analysis original dataset
write.csv(CCEP_amy_LDLPFC, paste0(rootpath, "TimeFrequency/CCEP_amy_ChlxTimes_L-DLPFC.csv"))
