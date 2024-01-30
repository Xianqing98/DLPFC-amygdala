library(tidyverse)
## Plot Individual Trials ##
# rootDir = "/Users/xianqing/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/"
rootDir = "/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/"
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_L-DLPFC_amy/1hp_detrend500")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_L-DLPFC_amy/2hp_detrend500")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_L-DLPFC_amy/2hp_detrend500_cutaneous")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_L-DLPFC_amy/2hp_detrend500_sgACC")

wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_R-DLPFC_amy/2hp_detrend500") ## TMS R DLPFC amy
wd <- paste0(rootDir,"iTEP_csv_R-DLPFC_amy/unfilter") 

wd <- paste0(rootDir,"allPatients_iEEGdata_csv/iTEP_csv_R-Parietal_amy/2hp_detrend500")## TMS R Parietal amy
# ccep
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_L-DLPFC_amy/2hp_nodetrend")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_L_DLPFC/sessions")

wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_R-DLPFC_amy") ## esTT R DLPFC amy
wd <- paste0(rootDir,"CCEP_csv_R-DLPFC_amy/unfilter") 

wd <- paste0(rootDir,"CCEP_csv_L-Parietal_amy") ## esTT L Parietal amy
wd <- paste0(rootDir,"CCEP_csv_L-Parietal_amy/unfilter") 

wd <- paste0(rootDir,"CCEP_csv_R-Parietal_amy") ## esTT R Parietal amy
wd <- paste0(rootDir,"CCEP_csv_R-Parietal_amy/unfilter") 

wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_L-IFG_amy/2hp_nodetrend")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_L-IFG_amy/2hp_nodetrend_Hipp")
wd <- paste0(rootDir,"allPatients_iEEGdata_csv/CCEP_csv_L_VLPFC/sessions")


wd <- paste0(rootDir,"CCEP_csv_R-IFG_amy") ## esTT R IFG amy
wd <- paste0(rootDir,"CCEP_csv_R-IFG_amy/unfilter") 

# wd <- paste0(rootDir,"iTEP_csv_L-DLPFC_amy/672_mask") ## TMS L DLPFC amy

setwd(wd)
TEP_amy <- read_csv(list.files(wd, pattern = '*.csv'))
names(TEP_amy) <- c("Patient","Condition","TrialNumber","Channel","Time","Amplitude")
names(TEP_amy) <- c("Patient","TrialNumber","Channel","Time","Amplitude")
names(TEP_amy) <- c("Patient","Session","TrialNumber","Channel","Time","Amplitude")
TEP_amy$Time <- TEP_amy$Time*1000
TEP_amy$Channel <- paste(TEP_amy$Patient, TEP_amy$Channel, sep = "_")
TEP_amy$Channel <- paste(TEP_amy$Session, TEP_amy$Channel, sep = "_")
TEP_amy <- filter(TEP_amy, Patient != 672 & 
                    Channel != "593_LFPx183" &
                    Channel != "593_LFPx184" &
                    Channel != "593_LFPx185" &
                    Channel != "430_LFPx14")

ChannelList <- TEP_amy %>% group_by(Channel) %>% 
  summarise(Patient = mean(Patient))
ChannelList <- ChannelList[[1]]

################################################################################
filter <- 2    # 0     1hp     2hp
detrend <- 0.5   # 0.    0.5s.   1s
################################################################################

if (filter == 0){
  filetype <- "_unfilter.png"
} else if (filter == 1){
  if (detrend == 0.5){filetype <- "_1hp_detrend500.png"} 
  else if (detrend == 1){filetype <- "_1hp_detrend1000.png"}
  else if (detrend == 0){filetype <- "_1hp_nodetrend.png"}
} else if (filter == 2){
  if (detrend == 0.5){filetype <- "_2hp_detrend500.png"} 
  else if (detrend == 1){filetype <- "_2hp_detrend1000.png"}
  else if (detrend == 0){filetype <- "_2hp_nodetrend.png"}
  
} 

# plot iTEP
for (chl in ChannelList){
  for (cond in c("sham","tms")){
    if (cond == "tms"){
      color <- "dodgerblue3"
    } else {
      color <- "gray40"
    }
    print(paste0("Plotting ",chl,"_noDetrend_",cond,filetype)) # mark: no detrend
    p <- TEP_amy %>% filter(Channel == chl & Condition == cond) %>% 
      ggplot(aes(Time, Amplitude))+
      stat_summary(col = color, geom="line")+
      annotate("rect", xmin=-10, xmax=25, ymin=-75, ymax=75, fill="gray", alpha = 0.7)+
      # annotate("rect", xmin=-10, xmax=25, ymin=-75, ymax=75, fill="gray")+
      geom_hline(yintercept=0, linetype="dotted", col = "darkgray")+
      xlim(-250,500)+
      # xlim(-250,450)+
      scale_y_continuous(limits=c(-100,100), breaks = c(-50,50))+
      ylab("Trial")+xlab("Time (ms)")+
      ggtitle(paste(chl,cond,sep=", "))+theme_classic()+
      theme(axis.text = element_text(color = "black", size = rel(1.2)))+
      theme(axis.title = element_text(color = "black", size = rel(1.4)))+
      facet_wrap(~TrialNumber)
      # facet_wrap(~TrialNumber,ncol=4,strip.position = "top")
    # print(p)
    ggsave(paste0(chl,"_",cond,filetype), 
           p, height = 10, width = 15, units = "in", dpi = 900)
    print(paste0(chl,"_",cond,filetype," saved"))
  }
}

# plot CCEP
for (chl in ChannelList){
# for (chl in ChannelList[c(11:26,37:54)]){
    print(paste0("Plotting ",chl,"_esTT.png"))
    p <- TEP_amy %>% filter(Channel == chl) %>% 
      ggplot(aes(Time, Amplitude))+
      stat_summary(col = "royalblue3", geom="line")+
      annotate("rect", xmin=-8, xmax=8, ymin=-75, ymax=75, fill="gray",alpha = 0.7)+
      geom_hline(yintercept=0, linetype="dotted", col = "darkgray")+
      xlim(-100,500)+
      scale_y_continuous(limits=c(-100,100), breaks = c(-50,50))+
      ylab("Trial")+xlab("Time (ms)")+
      ggtitle(paste(chl,"es-TT",sep=", "))+theme_classic()+
      theme(axis.text = element_text(color = "black", size = rel(1.2)))+
      theme(axis.title = element_text(color = "black", size = rel(1.4)))+
      facet_wrap(~TrialNumber)
    # print(p)
    ggsave(paste0(chl,"_esTT",filetype), 
           p, height = 10, width = 15, units = "in", dpi = 900)
    print(paste0(chl,"_esTT", filetype, " saved"))
}
