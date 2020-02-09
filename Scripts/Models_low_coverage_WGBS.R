library(dplyr)
library(ggplot2)

# Set folder (Data folder)
setwd("~/PBAT_low_coverage/Data")

###################### Plot margin of error theoretical model #####################

# Theoretical model
zstar = qnorm(.995) #This calculates z* at a 99.5% confidence interval
p = 0.25 #The approximate methylation level from this example
n = seq(50,30000,by=50) #Numbers of cytosine calls assayed
E = (zstar*sqrt((p*(1-p))/n))*100 #Margin of error calculation given above variables

n_points = length(n) # Number of points theoretical model
t_model = data.frame("Pred_Cs"=n,"ME"=E,"Group"=rep("Model",n_points))

# Plot theoretical model
plot_t_model = ggplot(t_model, aes(x = Pred_Cs, y = ME, color = Group)) +
  geom_point(alpha = 0.9, size = 2.5) +
  geom_line(data = t_model, aes(x = Pred_Cs, y = ME)) + 
  xlab("Number of Cs in CG context") + ylab("Margin of error (± %) 99.5% CI") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        axis.line = element_line(size = 0.5)) +
  geom_hline(yintercept=1,linetype="dashed", color="blue", size = 0.75) +
  geom_hline(yintercept=2.5,linetype="dashed", color="green", size = 0.75) +
  geom_hline(yintercept=5,linetype="dashed", color="red",size = 0.75) +
  geom_vline(xintercept=10000,linetype="dashed", color="black",size = 0.75) +
  ylim(0,25) + scale_color_manual(values=c('#000000'))

# 1000 x 750
plot_t_model

###################### Plot example data bovine blastocyst #####################

# Import cattle blastocyst data
cattle_data <- read.delim("Cattle_Blastocyst.txt")
cattle_data$rank = rank(-cattle_data$CsM_CpG,ties.method=c("random"))
cattle_data = cattle_data[order(cattle_data$rank),]

# Estimating margin of error for table
zstar = qnorm(.995) #Inputted here is a 99.5% confidence interval
p = cattle_data$CsM_CpG #Methylation proportion of bovine blastocysts in this example
p = p / 100
n=cattle_data$Total_CpG #This is the number of cytosine calls assayed
ME = (zstar*sqrt((p*(1-p))/n))*100

# Barplot
plot_example <- ggplot(data=cattle_data, aes(x=rank,y=CsM_CpG,fill=trt)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=CsM_CpG-ME, ymax=CsM_CpG+ME),width=0.5,position=position_dodge(0.1)) +
  xlab("Sample") + ylab("CG methylation (%)") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        axis.line = element_line(size = 0.5)) +
  ylim(0,100)

# 1000 x 750
plot_example

###################### Plot dispersion bootstraping model #####################

# Data 1 - CF
data_CF <- read.csv("CF_summary_bootstraping.txt",sep = "\t", stringsAsFactors = FALSE)
data_CF <- data_CF[!apply(is.na(data_CF) | data_CF == "", 1, all),]
data_CF$Total_Cs_CpG = data_CF$Total_CsM_CpG + data_CF$Total_CsU_CpG
data_CF$CsM_CpG = data_CF$Total_CsM_CpG * 100 / data_CF$Total_Cs_CpG
data_CF$Reads <- sapply((strsplit(as.character(data_CF$Sample),'_')),"[",2)

# Data 2 - TP
data_TP <- read.csv("TP_summary_bootstraping.txt",sep = "\t", stringsAsFactors = FALSE)
data_TP <- data_TP[!apply(is.na(data_TP) | data_TP == "", 1, all),]
data_TP$Total_Cs_CpG = data_TP$Total_CsM_CpG + data_TP$Total_CsU_CpG
data_TP$CsM_CpG = data_TP$Total_CsM_CpG * 100 / data_TP$Total_Cs_CpG
data_TP$Reads <- sapply((strsplit(as.character(data_TP$Sample),'_')),"[",2)

# Merge data
data_all <- rbind(cbind(data.frame("Group"=rep("Female",21000)),data_CF),
                  cbind(data.frame("Group"=rep("Male",21000)),data_TP))

# Plot dispersion
plot_dispersion = ggplot(data_all, aes(x = Total_Cs_CpG, y = CsM_CpG, color = Group)) +
  geom_point(alpha = 0.5) +
  xlab("Number of Cs in CG context") + ylab("CG methylation (%)") +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    axis.line = element_line(size = 0.5)
  ) + 
  scale_color_manual(values=c('#F8766D','#619CFF')) +
  scale_fill_manual(values=c('#F8766D','#619CFF')) +
  ylim(25,100) +
  geom_hline(yintercept=80.20,linetype="dashed", color="blue") + 
  geom_hline(yintercept=69.4,linetype="dashed", color="red")

# Plot (1000 x 750)
plot_dispersion

###################### Plot margin of error bootstraping model #####################

# Data CF
code_CF <- read.csv("CF_code.csv",stringsAsFactors = FALSE)
code_CF <- lapply(code_CF, as.character)
map = setNames(code_CF$Group,code_CF$Number)

data_CF$Pred_Cs <- as.integer(map[unlist(data_CF$Reads)])

data_CF_summary <- as.data.frame(data_CF %>% group_by(Pred_Cs) %>%
                                      summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                                quantile_99_75 = quantile(CsM_CpG, c(.9975)),
                                                quantile_00_25 = quantile(CsM_CpG, c(.0025)),
                                                dif = quantile_99_75 - quantile_00_25,
                                                ME = dif / 2))

data_CF_summary

# Data TP
code_TP <- read.csv("TP_code.csv",stringsAsFactors = FALSE)
code_TP <- lapply(code_TP, as.character)
map = setNames(code_TP$Group,code_TP$Number)

data_TP$Pred_Cs <- as.integer(map[unlist(data_TP$Reads)])

data_TP_summary <- as.data.frame(data_TP %>% group_by(Pred_Cs) %>%
                                   summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                             quantile_99_75 = quantile(CsM_CpG, c(.9975)),
                                             quantile_00_25 = quantile(CsM_CpG, c(.0025)),
                                             dif = quantile_99_75 - quantile_00_25,
                                             ME = dif / 2))

data_TP_summary

# Bind data
data_summary <- rbind(data_CF_summary,data_TP_summary)
data_summary$Group = c(rep("Female",21),rep("Male",21))
data_summary

# Theoretical model
zstar = qnorm(.995) #This calculates z* at a 99.5% confidence interval
p = 0.75 #The approximate methylation level from this example
n = seq(50,30000,by=50) #Numbers of cytosine calls assayed
E = (zstar*sqrt((p*(1-p))/n))*100 #Margin of error calculation given above variables

n_points = length(n) # Number of points theoretical model
t_model = data.frame("Pred_Cs"=n,"ME"=E,"Group"=rep("Model",n_points))

# Plot margin of error
plot_ME = ggplot(data_summary, aes(x = Pred_Cs, y = ME, color = Group)) +
  geom_point(alpha = 0.9, size = 4) +
  geom_line(data = t_model, aes(x = Pred_Cs, y = ME)) + 
  xlab("Number of Cs in CG context") + ylab("Margin of error (± %) 99.5% CI") +
  theme(axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    axis.line = element_line(size = 0.5)) +
  geom_hline(yintercept=1,linetype="dashed", color="blue", size = 0.75) +
  geom_hline(yintercept=2.5,linetype="dashed", color="green", size = 0.75) +
  geom_hline(yintercept=5,linetype="dashed", color="red",size = 0.75) +
  ylim(0,25) + scale_color_manual(values=c('#F8766D','#619CFF','#000000'))

# 1000 x 750
plot_ME

###################### Plot margin of error bootstraping model and first CG #####################

# Data CF first CG
data_CF_fCG <- read.csv("data_CF_fCG.csv")

data_CF_fCG$CsM_CpG <- data_CF_fCG$CpG_met*100/data_CF_fCG$Sum
data_CF_fCG$Sample <- as.character(data_CF_fCG$Sample)
data_CF_fCG$Code <- sapply(strsplit(data_CF_fCG$Sample,"_"),"[",5)

CG_means_CF <- as.data.frame(data_CF_fCG %>% group_by(Code) %>%
                               summarise(mean = mean(Sum)))

map <- setNames(CG_means_CF$mean,CG_means_CF$Code)

data_CF_fCG$Pred_Cs <- as.integer(map[unlist(data_CF_fCG$Code)])

data_CF_fCG_summary <- as.data.frame(data_CF_fCG %>% group_by(Pred_Cs) %>%
                                       summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                                 quantile_99_75 = quantile(CsM_CpG, c(.9975)),
                                                 quantile_00_25 = quantile(CsM_CpG, c(.0025)),
                                                 dif = quantile_99_75 - quantile_00_25,
                                                 ME = dif / 2))

# Data TP first CG
data_TP_fCG <- read.csv("data_TP_fCG.csv")

data_TP_fCG$CsM_CpG <- data_TP_fCG$CpG_met*100/data_TP_fCG$Sum
data_TP_fCG$Sample <- as.character(data_TP_fCG$Sample)
data_TP_fCG$Code <- sapply(strsplit(data_TP_fCG$Sample,"_"),"[",5)

CG_means_TP <- as.data.frame(data_TP_fCG %>% group_by(Code) %>%
                               summarise(mean = mean(Sum)))

map <- setNames(CG_means_TP$mean,CG_means_TP$Code)

data_TP_fCG$Pred_Cs <- as.integer(map[unlist(data_TP_fCG$Code)])

data_TP_fCG_summary <- as.data.frame(data_TP_fCG %>% group_by(Pred_Cs) %>%
                                       summarise(mean = mean(CsM_CpG, na.rm = TRUE),
                                                 quantile_99_75 = quantile(CsM_CpG, c(.9975)),
                                                 quantile_00_25 = quantile(CsM_CpG, c(.0025)),
                                                 dif = quantile_99_75 - quantile_00_25,
                                                 ME = dif / 2))

# Bind data
data_fCG_summary <- rbind(data_CF_fCG_summary,data_TP_fCG_summary)
data_fCG_summary$Group = c(rep("Female_first_CG",21),rep("Male_first_CG",21))

# Bind all models
data_complete <- rbind(data_summary,data_fCG_summary)

# Plot data complete
plot_complete = ggplot(data_complete, aes(x = Pred_Cs, y = ME, color = Group)) +
  geom_point(alpha = 0.9, size = 3) +
  geom_line(data = t_model, aes(x = Pred_Cs, y = ME)) + 
  xlab("Number of Cs in CG context") + ylab("Margin of error (± %) 99.5% CI") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        axis.line = element_line(size = 0.5)) +
  ylim(0,25) + scale_color_manual(values=c('#F8766D','#FFA500','#619CFF','#32CD32','#000000'))

# 1300 x 750
plot_complete