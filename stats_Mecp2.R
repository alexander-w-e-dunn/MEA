# main scriot for analysis and plotting
# 
# to add:
# compute cohens d save in table for each comparison like p values
# do non parametric mixed effects model
# 
# ########## firing activity ##########
# prelim data using old thresholding methods
# setwd("C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Manuscripts/Mecp2/2022_data/Preliminary results")
# new data using shuffle controls
setwd("C:/Users/alexd/OneDrive - University of Cambridge/Cam/PhD/Manuscripts/PhD_Thesis/Chapter_Mecp2_effect/Mecp2_jitter_10ms")
# setwd("C:/Users/alexd/OneDrive - University of Cambridge/Cam/PhD/Manuscripts/PhD_Thesis/Chapter_Mecp2_effect/Mecp2_jitter_10ms/prop10")
# create dataframe to store result of interaction
results_df_interaction <- data.frame(Metric = character(), 
                                     ATS_grp = numeric(), DF_grp = numeric(), p_grp = numeric(),
                                     ATS_div = numeric(), DF_div = numeric(), p_div = numeric(),
                                     ATS_int = numeric(), DF_int = numeric(), p_int = numeric(), stringsAsFactors = FALSE)

# load data (old data using propotion and abs thresholds)
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Manuscripts/Mecp2/2022_data/cwt_L-0.0627_Mecp2_spikes_Nbursts_STTC_data.csv',header = T)
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Manuscripts/Mecp2/2022_data/cSpikes_5SD_Mecp2_analysis_2022_Aug_50_ms_dt.csv',header = T)
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Manuscripts/Mecp2/2022_data/cSpikes_5SD_Mecp2_analysis_2022_Aug_20_ms_dt.csv',header = T)
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Manuscripts/Mecp2/2022_data/cSpikes_5SD_Mecp2_analysis_2022_Aug_10_ms_dt.csv',header = T)
# new data using shuffles
data <- read.csv(file ='C:/Users/alexd/OneDrive - University of Cambridge/Cam/PhD/Manuscripts/PhD_Thesis/Chapter_Mecp2_effect/Mecp2_jitter_10ms/cSpikes_-0p0627_Mecp2_analysis_2023_Feb_jitter_dt_10_ms_bursts_graph.csv',header = T)
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University of Cambridge/Cam/PhD/Manuscripts/PhD_Thesis/Chapter_Mecp2_effect/Mecp2_jitter_10ms/prop10/cSpikes_-0p0627_Mecp2_analysis_2023_Feb_jitter_dt_10_ms_graph_prop10.csv',header = T)
data <- data[,2:ncol(data)]

# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Project/Analysis_and_reports/First_year_report/Data/all_data_netwSpikes_BakkumBursts_RSbursts_Regularity.csv',header = T)
# head(data)
# firing rates, network bursts, connectivity and graph metrics
# data <- read.csv(file ='C:/Users/alexd/OneDrive - University Of Cambridge/Cam/PhD/Project/Analysis_and_reports/First_year_report/Data/All_results_fr_netwBursts_fc_graph.csv',header = T)
#data <- read.csv(file ='C:/Users/owner/OneDrive - University Of Cambridge/Cam/PhD/Project/Analysis_and_reports/First_year_report/Data/All_results_fr_netwBursts_fc_graph.csv',header = T)
# head(data)
data <- within(data, Genotype <- factor(Genotype, levels = c("WT","HE","KO")))
data.all <- data
data.all <- within(data.all, Genotype <- factor(Genotype, levels = c("WT","HE","KO")))

metrics_names <- colnames(data)[1:ncol(data)]
# number of columns before DVs or output metrics
num_IVs <- 4

# data.ko <- subset(data,data$Genotype == "WT" | data$Genotype == "KO",drop = TRUE)
# data.ko <- within(data.ko, Genotype <- factor(Genotype, levels = c("WT","KO")))
# data.he <- subset(data,data$Genotype == "WT" | data$Genotype == "HE")
# data.he <- within(data.he, Genotype <- factor(Genotype, levels = c("WT","HE")))

# data <- data.ko
data <- data.all

# add character class DIV factor
data$DIV.char <- as.character(data$DIV)
data <- within(data, DIV.char <- factor(DIV.char, levels = c("7","14","21","28","35")))

# plot
# install.packages("ggplot2")
# install.packages("dplyr")
# if(!require(devtools)) install.packages("devtools")
# install.packages("plotrix") 
# install.packages("car")
# install.packages("dunn.test")
# install.packages("ggsignif")
# install.packages("lifecycle")
# install.packages("plotrix")
# install.packages("ImpactEffectsize")
# install.packages("lsr")

library(car)
library(dunn.test)
library(plotrix) 
library(ggplot2)
library(dplyr)
library(ggsignif)
library(lifecycle)
library(plotrix)
library(nparLD)
library(ImpactEffectsize)
library(lsr)

### general settings
# X <- data$Age
X <- as.character(data$DIV)
gtcolours<-c( "blue3","magenta4","firebrick3")
gtlines2 <- c("solid","dotdash","dashed")
gtlines  <- c("solid","solid","solid")
xstring <- "Days in vitro"

for (metric in c((1+num_IVs):length(metrics_names))) {
# for (metric in c(76:length(metrics_names))) {
# for (metric in c(5:8)) {
    print(metrics_names[metric])
  # remove nan values
  data[is.na(data[,metric]),metric] <- 0
  data$dv <- data[,metric]
  Y <- data$dv
  ylims.min <- min(na.omit(data$dv))
  ylims.max <- max(na.omit(data$dv))
  titlestring <- metrics_names[metric]
  ystring <- metrics_names[metric]
  boxfilestring<- paste("MECP2_all_Box",metrics_names[metric],".png",sep = "_")
  boxfilestring2<- paste("MECP2_all_Box",metrics_names[metric],"_no_lines.png",sep = "_")
  linefilestring<- paste("MECP2_all_Lin",metrics_names[metric],".png",sep = "_")
  linefilestring2<- paste("MECP2_all_Lin_Cultures",metrics_names[metric],".png",sep = "_")
  
  
  ################# plotting
  gd <- data %>% 
    group_by(Genotype, DIV) %>% 
    summarise(dv = mean(na.omit(dv)))
  gd
  
  sd <- data %>%
    group_by(Genotype, DIV) %>%
    summarise(sdpergroup = sd(na.omit(dv)))
  sd
  
  nums <- data %>%
    group_by(Genotype, DIV) %>%
    summarise(numspergroup = n_distinct(Culture))
  nums
  
  gd$sd <- sd$sdpergroup
  gd$ns <- nums$numspergroup
  gd$sem <- gd$sd / (sqrt(gd$ns))
  
  gd$dv <- round(gd$dv,digits = 3)
  gd$sd <- round(gd$sd,digits = 3)
  gd$sem <- round(gd$sem,digits = 3)
  
  gd
  
  subdata <- subset(data,!is.na(data$dv))
  subdata.all <- subdata
  aggregate(dv ~ DIV:Genotype, data = subdata, 
            FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x), n = length(x)))
  
  # output m and sem in table and save
  # capture.output(data.frame(gd), file = paste("MECP2_all_M_SEM",metrics_names[metric],".txt",sep = "_"))
  # write.csv(data.frame(gd), file = paste("MECP2_all_M_SEM",metrics_names[metric],".txt",sep = "_"))
  write.csv(data.frame(gd), file = paste("MECP2_all_M_SEM",metrics_names[metric],".csv",sep = "_"),row.names = FALSE)
  # output corrected p values
  
  # line plot
  # png(file=linefilestring,width=1200,height=800,res=300)
  
  ggplot(data, aes(x = DIV, y = dv, color = Genotype,linetype = Genotype)) +
    # geom_point(aes(group = Culture), alpha = .2,size=0.5,position=position_dodge(1)) +
    # geom_line(aes(group = Culture), alpha = .2,size=0.5,position=position_dodge(1)) +
    # geom_line(aes(linetype = Genotype),data = gd, alpha = 1, size = 0.7) +
    geom_line(aes(group = Genotype), data = gd, alpha = 1, size = .75,position=position_dodge(3.5)) +
    scale_linetype_manual(breaks = c("WT","HE", "KO"), values=gtlines) +
    geom_errorbar(data = gd, aes(ymin=dv-sem, ymax=dv+sem,group = Genotype), width=1.1,
                  position=position_dodge(3.5),size = 0.75,linetype = "solid") +
    # geom_errorbar(df2,aes(ymin=array_fracInBursts-sd, ymax=array_fracInBursts+sd,group = Genotype), width=.2,
    #               position=position_dodge(0.05)) 
    # 
    geom_point(aes(group = Genotype,pch = Genotype), data = gd, alpha = 1, size = 1.5,position=position_dodge(3.5)) +
    scale_x_continuous(name = xstring, breaks = c(7,14,21,28,35),labels = c("7","14","21","28","35")) +
    scale_color_manual(name = "Genotype",breaks = c("WT","HE","KO"), values=gtcolours) +
    # ylim(ylims.min,ylims.max) +
    theme_classic() +
    theme(text = element_text(size=14)) +
    labs(
      title = titlestring,
      x = xstring,
      y = ystring,
      color = NULL
    )
  
  ggsave(
    filename = linefilestring,
    plot = last_plot(),
    device = "png",
    scale = 8,
    width = 1.5,
    height = 1,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
    )   
  
  # dev.off()
  
  #  boxplot
  
  # png(file=boxfilestring,width=1200,height=800,res=300)
  
  ggplot(data, aes(DIV.char, dv , colour = Genotype), group = Genotype) + 
    # could add the ablines from linear model given slope and intercept for each group; could be done manually with geomabline or with geomsmooth
    # geom_abline(aes(group = Genotype), alpha = .3,size=0.5,linetype = "solid") +
    # geom_smooth(se = FALSE, method = lm) +
    # geom_smooth(aes(colour=Genotype)method = "nls", formula = y ~ a * x + b, se = F) +
    # geom_line(aes(group = Culture), alpha = .2,size=0.25,position=position_dodge(width = .75)) +
    geom_path(aes(group = Genotype), position = position_dodge(width = .75), size = .3, alpha = .1) +
    geom_boxplot(mapping = aes(x = DIV.char, y = Y,fill = Genotype),alpha = .5,outlier.size = 0,size=.3,outlier.shape = NA) +
    scale_x_discrete(name = xstring, breaks = c("7","14","21","28","35"),labels = c("7","14","21","28","35")) +
    scale_color_manual(breaks = c("WT", "HE", "KO"),
                       values=c("black","black","black"),labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) + 
    scale_fill_manual(breaks = c("WT", "HE","KO"),values=gtcolours,labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) +
    scale_shape_manual(breaks = c("WT", "HE","KO"),values=c(0,7,4),labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) +
    # scale_shape_manual(breaks = c("WT", "HE","KO"),values=c(16,17,4)) +
    geom_point(aes(pch = Genotype), position = position_jitterdodge(jitter.width=0.00,dodge.width=.75),size = 1.2) +
    # scale_size_manual(breaks = c("WT", "HE","KO"),values=c(8,1,1)) +
    ylim(ylims.min,ylims.max) +
    theme_classic() +
    theme(text = element_text(size=14),axis.text.x = element_text(color = "black"),
          legend.text.align = 0, axis.text.y = element_text(color = "black")) +
    
    labs(
      title = titlestring,
      y = ystring
    ) 
  
  ggsave(
    filename = boxfilestring,
    plot = last_plot(),
    device = "png",
    scale = 8,
    width = 1.5,
    height = 1,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )   
  
  ggplot(data, aes(DIV.char, dv , colour = Genotype), group = Genotype) + 
    # could add the ablines from linear model given slope and intercept for each group; could be done manually with geomabline or with geomsmooth
    # geom_abline(aes(group = Genotype), alpha = .3,size=0.5,linetype = "solid") +
    # geom_smooth(se = FALSE, method = lm) +
    # geom_smooth(aes(colour=Genotype)method = "nls", formula = y ~ a * x + b, se = F) +
    # geom_line(aes(group = Culture), alpha = .2,size=0.25,position=position_dodge(width = .75)) +
    # geom_path(aes(group = Genotype), position = position_dodge(width = .75), size = .3, alpha = .1) +
    geom_boxplot(mapping = aes(x = DIV.char, y = Y,fill = Genotype),alpha = .5,outlier.size = 0,size=.3,outlier.shape = NA) +
    scale_x_discrete(name = xstring, breaks = c("7","14","21","28","35"),labels = c("7","14","21","28","35")) +
    scale_color_manual(breaks = c("WT", "HE", "KO"),
                       values=c("black","black","black"),labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) + 
    scale_fill_manual(breaks = c("WT", "HE","KO"),values=gtcolours,labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) +
    scale_shape_manual(breaks = c("WT", "HE","KO"),values=c(0,7,4),labels=c("Wildtype", expression(paste(italic("Mecp2")," Het")),expression(paste(italic("Mecp2")," KO")))) +
    # scale_shape_manual(breaks = c("WT", "HE","KO"),values=c(16,17,4)) +
    geom_point(aes(pch = Genotype), position = position_jitterdodge(jitter.width=0.25,dodge.width=.75),size = 1.2) +
    # scale_size_manual(breaks = c("WT", "HE","KO"),values=c(8,1,1)) +
    ylim(ylims.min,ylims.max) +
    theme_classic() +
    theme(text = element_text(size=14),axis.text.x = element_text(color = "black"),
          legend.text.align = 0, axis.text.y = element_text(color = "black")) +
    
    labs(
      title = titlestring,
      y = ystring
    ) 
  
  ggsave(
    filename = boxfilestring2,
    plot = last_plot(),
    device = "png",
    scale = 8,
    width = 1.5,
    height = 1,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )   
  
  # dev.off()
  
  # line plot with all cultures
  # png(file=linefilestring2,width=1200,height=800,res=300)
  
  ggplot(data, aes(x = DIV, y = dv, color = Genotype,linetype = Genotype)) +
    geom_line(aes(group = Culture), alpha = .2,size=0.5,position=position_dodge(0.1)) +
    # geom_line(aes(linetype = Genotype),data = gd, alpha = 1, size = 0.7) +
    geom_line(aes(group = Genotype), data = gd, alpha = 1, size = 0.75,position=position_dodge(3.5)) +
    scale_linetype_manual(breaks = c("WT","HE", "KO"), values=gtlines2) +
    geom_errorbar(data = gd, aes(ymin=dv-sem, ymax=dv+sem,group = Genotype), width=1.1,
                  position=position_dodge(3.5),size = 0.75,linetype = "solid") +
    # geom_errorbar(df2,aes(ymin=array_fracInBursts-sd, ymax=array_fracInBursts+sd,group = Genotype), width=.2,
    #               position=position_dodge(0.05)) 
    # 
    scale_x_continuous(name = xstring, breaks = c(7,14,21,28,35),labels = c("7","14","21","28","35")) +
    scale_color_manual(name = "Genotype",breaks = c("WT","HE","KO"), values=gtcolours) +
    ylim(ylims.min,ylims.max) +
    theme_classic() +
    theme(text = element_text(size=14)) +
    labs(
      title = titlestring,
      x = xstring,
      y = ystring,
      color = NULL
    )
  
  ggsave(
    filename = linefilestring2,
    plot = last_plot(),
    device = "png",
    scale = 8,
    width = 1.5,
    height = 1,
    units = "cm",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
    )  
  
  # dev.off()
  
  ######### do friedman rep measures test on DIV and k-wallis test on genotype if anova doesnt meet assumptions
  
  aggregate(dv~DIV:Genotype,data = subdata,summary)
  aggregate(dv ~ DIV:Genotype, data = subdata, 
            FUN = function(x) c(mean = mean(x), se = std.error(x), sd = sd(x), n = length(x)))
  gd
  d07 <- subset(subdata,subdata$DIV == 7)
  d14 <- subset(subdata,subdata$DIV == 14)
  d21 <- subset(subdata,subdata$DIV == 21)
  d28 <- subset(subdata,subdata$DIV == 28)
  d35 <- subset(subdata,subdata$DIV == 35)
  
  d07all <- subset(subdata.all,subdata.all$DIV == 7)
  d14all <- subset(subdata.all,subdata.all$DIV == 14)
  d21all <- subset(subdata.all,subdata.all$DIV == 21)
  d28all <- subset(subdata.all,subdata.all$DIV == 28)
  d35all <- subset(subdata.all,subdata.all$DIV == 35)
  
  kw07 <- kruskal.test(dv~Genotype, data=d07all)
  kw14 <- kruskal.test(dv~Genotype, data=d14all)
  kw21 <- kruskal.test(dv~Genotype, data=d21all)
  kw28 <- kruskal.test(dv~Genotype, data=d28all)
  kw35 <- kruskal.test(dv~Genotype, data=d35all)
  kw_ps <- c(kw07$p.value, kw14$p.value, kw21$p.value, kw28$p.value, kw35$p.value)
  adj.kw.ps<-p.adjust(kw_ps,method = "fdr")
  
  mw07 <- wilcox.test(subset(d07$dv,d14$Genotype=="WT"), subset(d07$dv,d14$Genotype=="HE"), paired = FALSE, alternative = "two.sided")
  mw14 <- wilcox.test(subset(d14$dv,d14$Genotype=="WT"), subset(d14$dv,d14$Genotype=="HE"), paired = FALSE, alternative = "two.sided")
  mw21 <- wilcox.test(subset(d21$dv,d14$Genotype=="WT"), subset(d21$dv,d14$Genotype=="HE"), paired = FALSE, alternative = "two.sided")
  mw28 <- wilcox.test(subset(d28$dv,d14$Genotype=="WT"), subset(d28$dv,d14$Genotype=="HE"), paired = FALSE, alternative = "two.sided")
  mw35 <- wilcox.test(subset(d35$dv,d14$Genotype=="WT"), subset(d35$dv,d14$Genotype=="HE"), paired = FALSE, alternative = "two.sided")
  mw_ps <- c(mw07$p.value, mw14$p.value, mw21$p.value, mw28$p.value, mw35$p.value) # all mann whitney p values
  adj.mw.ps<-p.adjust(mw_ps,method = "fdr")
  
  mw07 <- wilcox.test(subset(d07$dv,d14$Genotype=="WT"), subset(d07$dv,d14$Genotype=="KO"), paired = FALSE, alternative = "two.sided")
  mw14 <- wilcox.test(subset(d14$dv,d14$Genotype=="WT"), subset(d14$dv,d14$Genotype=="KO"), paired = FALSE, alternative = "two.sided")
  mw21 <- wilcox.test(subset(d21$dv,d14$Genotype=="WT"), subset(d21$dv,d14$Genotype=="KO"), paired = FALSE, alternative = "two.sided")
  mw28 <- wilcox.test(subset(d28$dv,d14$Genotype=="WT"), subset(d28$dv,d14$Genotype=="KO"), paired = FALSE, alternative = "two.sided")
  mw35 <- wilcox.test(subset(d35$dv,d14$Genotype=="WT"), subset(d35$dv,d14$Genotype=="KO"), paired = FALSE, alternative = "two.sided")
  mw_ps <- c(mw07$p.value, mw14$p.value, mw21$p.value, mw28$p.value, mw35$p.value) # all mann whitney p values
  adj.mw.ps<-p.adjust(mw_ps,method = "fdr")
  
  d07uns <- unstack(d07all[c(length(colnames(subdata.all)),4)])
  d14uns <- unstack(d14all[c(length(colnames(subdata.all)),4)])
  d21uns <- unstack(d21all[c(length(colnames(subdata.all)),4)])
  d28uns <- unstack(d28all[c(length(colnames(subdata.all)),4)])
  d35uns <- unstack(d35all[c(length(colnames(subdata.all)),4)])
  
  dt07<-dunn.test(d07uns)
  dt14<-dunn.test(d14uns)
  dt21<-dunn.test(d21uns)
  dt28<-dunn.test(d28uns)
  dt35<-dunn.test(d35uns)
  
  gd.all <- gd
  
  subset(gd.all,gd.all$DIV == 07)
  subset(gd.all,gd.all$DIV == 14)
  subset(gd.all,gd.all$DIV == 21)
  subset(gd.all,gd.all$DIV == 28)
  subset(gd.all,gd.all$DIV == 35)
  
  c(round(mw07$statistic,3),round(adj.mw.ps[1],3),round(mw_ps[1],3))
  c(round(mw14$statistic,3),round(adj.mw.ps[1],3),round(mw_ps[1],3))
  c(round(mw21$statistic,3),round(adj.mw.ps[2],3),round(mw_ps[2],3))
  c(round(mw28$statistic,3),round(adj.mw.ps[3],3),round(mw_ps[3],3))
  c(round(mw35$statistic,3),round(adj.mw.ps[4],3),round(mw_ps[4],3))
  
  c(round(kw07$statistic,3),round(adj.kw.ps[1],3),round(kw_ps[1],3))
  c(round(kw14$statistic,3),round(adj.kw.ps[1],3),round(kw_ps[1],3))
  c(round(kw21$statistic,3),round(adj.kw.ps[2],3),round(kw_ps[2],3))
  c(round(kw28$statistic,3),round(adj.kw.ps[3],3),round(kw_ps[3],3))
  c(round(kw35$statistic,3),round(adj.kw.ps[4],3),round(kw_ps[4],3))
  
  dt14$comparisons
  round(dt07$P.adjusted,3)
  round(dt14$P.adjusted,3)
  round(dt21$P.adjusted,3)
  round(dt28$P.adjusted,3)
  round(dt35$P.adjusted,3)
  
  # Impact comparisons
  # DIV7 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D07.impact <- 
    c(Impact(subset(d07$dv,d07$Genotype=="WT"|d07$Genotype=="HE"),
             subset(d07$Genotype,d07$Genotype=="WT"|d07$Genotype=="HE"),PlotIt = FALSE)$Impact,
      Impact(subset(d07$dv,d07$Genotype=="WT"|d07$Genotype=="KO"),
             subset(d07$Genotype,d07$Genotype=="WT"|d07$Genotype=="KO"),PlotIt = FALSE)$Impact,
      Impact(subset(d07$dv,d07$Genotype=="HE"|d07$Genotype=="KO"),
             subset(d07$Genotype,d07$Genotype=="HE"|d07$Genotype=="KO"),PlotIt = FALSE)$Impact
  )
  
  # DIV14
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D14.impact <- 
    c(Impact(subset(d14$dv,d14$Genotype=="WT"|d14$Genotype=="HE"),
             subset(d14$Genotype,d14$Genotype=="WT"|d14$Genotype=="HE"),PlotIt = FALSE)$Impact,
      Impact(subset(d14$dv,d14$Genotype=="WT"|d14$Genotype=="KO"),
             subset(d14$Genotype,d14$Genotype=="WT"|d14$Genotype=="KO"),PlotIt = FALSE)$Impact,
      Impact(subset(d14$dv,d14$Genotype=="HE"|d14$Genotype=="KO"),
             subset(d14$Genotype,d14$Genotype=="HE"|d14$Genotype=="KO"),PlotIt = FALSE)$Impact
    )
  
  # DIV21
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D21.impact <- 
    c(Impact(subset(d21$dv,d21$Genotype=="WT"|d21$Genotype=="HE"),
             subset(d21$Genotype,d21$Genotype=="WT"|d21$Genotype=="HE"),PlotIt = FALSE)$Impact,
      Impact(subset(d21$dv,d21$Genotype=="WT"|d21$Genotype=="KO"),
             subset(d21$Genotype,d21$Genotype=="WT"|d21$Genotype=="KO"),PlotIt = FALSE)$Impact,
      Impact(subset(d21$dv,d21$Genotype=="HE"|d21$Genotype=="KO"),
             subset(d21$Genotype,d21$Genotype=="HE"|d21$Genotype=="KO"),PlotIt = FALSE)$Impact
    )
  
  # DIV28
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D28.impact <- 
    c(Impact(subset(d28$dv,d28$Genotype=="WT"|d28$Genotype=="HE"),
             subset(d28$Genotype,d28$Genotype=="WT"|d28$Genotype=="HE"),PlotIt = FALSE)$Impact,
      Impact(subset(d28$dv,d28$Genotype=="WT"|d28$Genotype=="KO"),
             subset(d28$Genotype,d28$Genotype=="WT"|d28$Genotype=="KO"),PlotIt = FALSE)$Impact,
      Impact(subset(d28$dv,d28$Genotype=="HE"|d28$Genotype=="KO"),
             subset(d28$Genotype,d28$Genotype=="HE"|d28$Genotype=="KO"),PlotIt = FALSE)$Impact
    )
  
  # DIV35
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D35.impact <- 
    c(Impact(subset(d35$dv,d35$Genotype=="WT"|d35$Genotype=="HE"),
             subset(d35$Genotype,d35$Genotype=="WT"|d35$Genotype=="HE"),PlotIt = FALSE)$Impact,
      Impact(subset(d35$dv,d35$Genotype=="WT"|d35$Genotype=="KO"),
             subset(d35$Genotype,d35$Genotype=="WT"|d35$Genotype=="KO"),PlotIt = FALSE)$Impact,
      Impact(subset(d35$dv,d35$Genotype=="HE"|d35$Genotype=="KO"),
             subset(d35$Genotype,d35$Genotype=="HE"|d35$Genotype=="KO"),PlotIt = FALSE)$Impact
    )
  
  # COHENS D CALCULATIONS
  # DIV07 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D07.cohen <- 
    c(cohensD(subset(d07$dv,d07$Genotype=="WT"),subset(d07$dv,d07$Genotype=="HE")),
      cohensD(subset(d07$dv,d07$Genotype=="WT"),subset(d07$dv,d07$Genotype=="KO")),
      cohensD(subset(d07$dv,d07$Genotype=="HE"),subset(d07$dv,d07$Genotype=="KO"))
    )
  # DIV14 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D14.cohen <- 
    c(cohensD(subset(d14$dv,d14$Genotype=="WT"),subset(d14$dv,d14$Genotype=="HE")),
      cohensD(subset(d14$dv,d14$Genotype=="WT"),subset(d14$dv,d14$Genotype=="KO")),
      cohensD(subset(d14$dv,d14$Genotype=="HE"),subset(d14$dv,d14$Genotype=="KO"))
    )
  # DIV21 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D21.cohen <- 
    c(cohensD(subset(d21$dv,d21$Genotype=="WT"),subset(d21$dv,d21$Genotype=="HE")),
      cohensD(subset(d21$dv,d21$Genotype=="WT"),subset(d21$dv,d21$Genotype=="KO")),
      cohensD(subset(d21$dv,d21$Genotype=="HE"),subset(d21$dv,d21$Genotype=="KO"))
    )
  # DIV28 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D28.cohen <- 
    c(cohensD(subset(d28$dv,d28$Genotype=="WT"),subset(d28$dv,d28$Genotype=="HE")),
      cohensD(subset(d28$dv,d28$Genotype=="WT"),subset(d28$dv,d28$Genotype=="KO")),
      cohensD(subset(d28$dv,d28$Genotype=="HE"),subset(d28$dv,d28$Genotype=="KO"))
    )
  # DIV35 
  # WT vs HE 
  # WT vs KO
  # HE vs KO
  D35.cohen <- 
    c(cohensD(subset(d35$dv,d35$Genotype=="WT"),subset(d35$dv,d35$Genotype=="HE")),
      cohensD(subset(d35$dv,d35$Genotype=="WT"),subset(d35$dv,d35$Genotype=="KO")),
      cohensD(subset(d35$dv,d35$Genotype=="HE"),subset(d35$dv,d35$Genotype=="KO"))
    )
  
  
  df<-data.frame("K","p","p_adj_fdr")
  names(df)<-c("KW statistic","KW p value","KW FDR-adjusted p")
  df[2,]<- matrix(c(round(kw07$statistic,3),round(kw_ps[1],3),round(adj.kw.ps[1],3)),1,3)
  df[3,]<- matrix(c(round(kw14$statistic,3),round(kw_ps[2],3),round(adj.kw.ps[2],3)),1,3)
  df[4,]<- matrix(c(round(kw21$statistic,3),round(kw_ps[3],3),round(adj.kw.ps[3],3)),1,3)
  df[5,]<- matrix(c(round(kw28$statistic,3),round(kw_ps[4],3),round(adj.kw.ps[4],3)),1,3)
  df[6,]<- matrix(c(round(kw35$statistic,3),round(kw_ps[5],3),round(adj.kw.ps[5],3)),1,3)
  df[7,]<- data.frame("dunn test p","dunn test p","dunn test p")
  # df[8,]<- data.frame(t(dt14$comparisons))
  df[8,]<- data.frame(t(c("WT vs Mecp2-het","WT vs Mecp2_KO","Mecp2-Het vs Mecp2-KO")))
  df[9,]<- round(dt07$P.adjusted,3)
  df[10,]<- round(dt14$P.adjusted,3)
  df[11,]<- round(dt21$P.adjusted,3)
  df[12,]<- round(dt28$P.adjusted,3)
  df[13,]<- round(dt35$P.adjusted,3)
  
  # add cohen's d (and ImpactEffectSize due to non normal data Lötsch and Ultsch 2020)
  df[14,]<- data.frame(""," "," ")
  df[15,]<- round(D07.cohen,3)
  df[16,]<- round(D14.cohen,3)
  df[17,]<- round(D21.cohen,3)
  df[18,]<- round(D28.cohen,3)
  df[19,]<- round(D35.cohen,3)
  
  df[20,]<- data.frame(" "," "," ")
  df[21,]<- round(D07.impact,3)
  df[22,]<- round(D14.impact,3)
  df[23,]<- round(D21.impact,3)
  df[24,]<- round(D28.impact,3)
  df[25,]<- round(D35.impact,3)
  
  rnames <- c("DIV07","DIV14","DIV21","DIV28","DIV35",
              "Dunn's test","DIV07","DIV14","DIV21","DIV28","DIV35",
              "Cohen's d","DIV07","DIV14","DIV21","DIV28","DIV35",
              "Impact","DIV07","DIV14","DIV21","DIV28","DIV35")
  
  write.csv(df[c(2:6, 8:25),],file = paste("MECP2_stats_tests_",metrics_names[metric],".csv",sep = "_"),
            row.names = rnames)
  
  ats.int<-nparLD(dv ~ Genotype * DIV, data = data, subject = "Culture", description = FALSE)
  summary(ats.int)
  sm.int<-summary(ats.int)
  results_df_interaction <- rbind(results_df_interaction, data.frame(Metric = metrics_names[metric], 
                                                                     ATS_grp = sm.int$ANOVA.test[1,1], 
                                                                     DF_grp = sm.int$ANOVA.test[1,2], 
                                                                     p_grp = sm.int$ANOVA.test[1,3],
                                                                     ATS_div = sm.int$ANOVA.test[2,1], 
                                                                     DF_div = sm.int$ANOVA.test[2,2], 
                                                                     p_div = sm.int$ANOVA.test[2,3],
                                                                     ATS_int = sm.int$ANOVA.test[3,1], 
                                                                     DF_int = sm.int$ANOVA.test[3,2], 
                                                                     p_int = sm.int$ANOVA.test[3,3]))
  ### save p values and outputs like MW-U and KW stat 
  ### add friedman test between time points and add stars to plot
  
  
}

write.csv(results_df_interaction, file = "Mecp2_dt10ms_interaction_results_connections_dt10ms_edgelength.csv", row.names = FALSE)


