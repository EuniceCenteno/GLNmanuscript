##----        R-Script for piglet ileal microbiota analysis 
#                      Alpha diversity statistical analysis
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020

### 
library(stats)
library(ggplot2)
library(fdrtool)
library(ggpubr)
library(tidyverse)
library(car) ##for type III SS 
library(dplyr)
library(emmeans)
library(sjstats)
library(ggfortify)
library(summarytools)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)

setwd("~/Desktop/eunice/Glutamine/No.D0/metadata/")
#alpha diversity dataset
alpha_div <- read.table(file = "N0.list.stability.opti_mcc.0.03.pick.0.03.subsample.groups.summary.txt", sep = "\t", header = T)
#dataset
design <- read.table(file = "No0.JayJohnson2.design.txt", sep = "\t", header = TRUE)
str(design)
design <- design[,-c(1,6:7)] #remove unnecessary data
alpha_div_merge <- merge(alpha_div, design, by.x = "Group", by.y = "Group")
str(alpha_div_merge)

#Day 13
alpha_d13 <- subset(alpha_div_merge, day == 13 )
str(alpha_d13)
alpha_d13$Group <- factor(alpha_d13$Group)
alpha_d13$day <- factor(alpha_d13$day)
alpha_d13$season <- factor(alpha_d13$season)
alpha_d13$treatment <- factor(alpha_d13$treatment)
levels(alpha_d13$season) <- list("Summer"="Summer", "Spring"="Spring")
levels(alpha_d13$treatment) <- list("AB"="A", "GLN"="GLN", "NA"="NA_C")

#Day 33
alpha_d33 <- subset(alpha_div_merge, day == 33)
str(alpha_d33)
alpha_d33$Group <- factor(alpha_d33$Group)
alpha_d33$day <- factor(alpha_d33$day)
alpha_d33$season <- factor(alpha_d33$season)
alpha_d33$treatment <- factor(alpha_d33$treatment)
levels(alpha_d33$season) <- list("Summer"="Summer", "Spring"="Spring")
levels(alpha_d33$treatment) <- list("AB"="A", "GLN"="GLN", "NA"="NA_C")


### Checking is the data is normally distributed 
shapiro.test(alpha_d13$chao) 
shapiro.test(alpha_d13$shannon)

shapiro.test(alpha_d33$chao) #not normal
shapiro.test(alpha_d33$shannon)


# ----Defining GLM models for normal distributed data
#----Day 13: Chao and Shannon
anov_D13chao <- lm(chao ~ treatment + season + treatment:season, data = alpha_d13) 
anov_D13sha <- lm(shannon ~ treatment + season + treatment:season, data = alpha_d13) 

#----Day 33: Shannon
anov_D33chao <- lm(chao ~ treatment + season + treatment:season, data = alpha_d33) 
anov_D33sha <- lm(shannon ~ treatment + season + treatment:season, data = alpha_d33) 

##Kruskal-wallis for non-normal distributed data
kruskal.test(chao ~ treatment, data=alpha_d33)
kruskal.test(chao ~ season, data=alpha_d33)

interTS<-interaction(alpha_d33$treatment, alpha_d33$season)
kruskal.test(chao ~ interTS, data=alpha_d33)
library(rstatix)
alpha_d33 %>% kruskal_effsize(chao ~ interTS)


##---------------- ANOVA ANALYSIS------------- # 
str(alpha_div_merge)
#chao and shannon are number
alpha_div_merge %>% tally()
alpha_div_merge %>% count(day)
#d13    56
#d33    51

## ANOVA TYPE III SUM OF SQUARE
#you used the lm model
options(contrasts=c("contr.sum", "contr.poly"))

###----------------------------D 13 ---------------------# 

# DEFAULT ANOVA
#Create the model
#--------------- CHAO -------------
#1. OTUS COUNTS

#Anova type III SS
Anov3_D13chao <- Anova(anov_D13chao, type = 3)
Anov3_D13chao

#### calculating effect size, partial-eta square 
# Using as reference the tutorial "ANOVA tables in R' by Rose Hartman, 2018
# http://www.understandingdata.net/2017/05/11/anova-tables-in-r/

#It’s the sums of squares for each effect divided by the error SS
sschaoD13 <- Anova(anov_D13chao, type = 3)
sschaoD13$pes <- c(sschaoD13$'Sum Sq'[-nrow(sschaoD13)], NA)/(sschaoD13$'Sum Sq' + sschaoD13$'Sum Sq'[nrow(sschaoD13)])
sschaoD13


#------------------SHANNON -----------
# 2. Shannon (Alpha diversity richness and evenness)
#you used the lm model
options(contrasts=c("contr.sum", "contr.poly"))

#ANOVA type III SS
Anov3_D13sha <- Anova(anov_D13sha, type = 3)
Anov3_D13sha

ssshaD13 <- Anova(anov_D13sha, type = 3)
ssshaD13$pes <- c(ssshaD13$'Sum Sq'[-nrow(ssshaD13)], NA)/(ssshaD13$'Sum Sq' + ssshaD13$'Sum Sq'[nrow(ssshaD13)])
ssshaD13

##----------------------------- DAY 33-------------------#

#Anova type III SS
#------------------SHANNON -----------
# 2. Shannon (Alpha diversity richness and evenness)
#you used the lm model
options(contrasts=c("contr.sum", "contr.poly"))

Anov3_D33sha <- Anova(anov_D33sha, type = 3)
Anov3_D33sha

sshaD33 <- Anova(anov_D33sha, type = 3)
sshaD33$pes <- c(sshaD33$'Sum Sq'[-nrow(sshaD33)], NA)/(sshaD33$'Sum Sq' + sshaD33$'Sum Sq'[nrow(sshaD33)])
sshaD33

#-----CALCULATION OF ERROR BARS
#use the model you created to run the ANOVAs
library(emmeans)
anov_D13chao 
anov_D13sha 

#----Day 33: Chao and Shannon
anov_D33chao 
anov_D33sha 


#Day 13 ~ chao
trt13_means <- emmeans(anov_D13chao, ~ treatment)
trt13_means_df <- data.frame(trt13_means)
trt13_means_df
day <- "Day 13"
trt13_means_df$day <- day

season13_means <- emmeans(anov_D13chao, ~ season)
season13_means_df <- data.frame(season13_means)
season13_means_df
day <- "Day 13"
season13_means_df$day <- day

#Day 13 ~ shannon
trt13_meanSha <- emmeans(anov_D13sha, ~ treatment)
trt13_means_dfSha <- data.frame(trt13_meanSha)
trt13_means_dfSha
day <- "Day 13"
trt13_means_dfSha$day <- day

sea13_meanSha <- emmeans(anov_D13sha, ~ season)
sea13_means_dfSha <- data.frame(sea13_meanSha)
sea13_means_dfSha
day <- "Day 13"
sea13_means_dfSha$day <- day

#Day 33
#Chao
trt33_means <- emmeans(anov_D33chao, ~ treatment)
trt33_means_df <- data.frame(trt33_means)
trt33_means_df
day <- "Day 33"
trt33_means_df$day <- day

season33_means <- emmeans(anov_D33chao, ~ season)
season33_means_df <- data.frame(season33_means)
season33_means_df
day <- "Day 33"
season33_means_df$day <- day

#shannon
trt33_meanSha <- emmeans(anov_D33sha, ~ treatment)
trt33_means_dfSha <- data.frame(trt33_meanSha)
trt33_means_dfSha
day <- "Day 33"
trt33_means_dfSha$day <- day

sea33_meanSha<-emmeans(anov_D33sha, ~ season)
sea33_means_dfSha <- data.frame(sea33_meanSha)
sea33_means_dfSha
day1 <- "Day 33"
sea33_means_dfSha$day <- day1


#--------FIGURES CHAO AND SHANNON------------#

my_colors2 <- c('black',"gray") 

##Chao d13
library(extrafont)
font_import()
loadfonts(device = "win")
fonts()

my_colors <- c("black","gray", "white") 
a  <- ggplot(data = trt13_means_df, aes(x = treatment, y = emmean, fill = treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Treatment") + ylab("Chao") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  theme(text=element_text(family="Times New Roman"))

w <- ggplot(data = season13_means_df, aes(x = season, y = emmean, fill = season)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Season") + ylab("Chao") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color="black", size=10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))

##Chao d33
c <- ggplot(data = trt33_means_df, aes(x = treatment, y = emmean, fill = treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Treatment") + ylab("Chao") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))


y <- ggplot(data = season33_means_df, aes(x = season, y = emmean, fill = season)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Season") + ylab("Chao") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank())+
  theme(text=element_text(family="Times New Roman"))

#--------------------------- shannon 13
b <- ggplot(data = trt13_means_dfSha, aes(x = treatment, y = emmean, fill = treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Treatment") + ylab("Shannon") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))


x <- ggplot(data = sea13_means_dfSha, aes(x = season, y = emmean, fill = season)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Season") + ylab("Shannon") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, , face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10))+
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))

#day33
d <- ggplot(data = trt33_means_dfSha, aes(x = treatment, y = emmean, fill = treatment)) + 
  geom_bar(stat = "identity", alpha=0.8, colour="black") +
  theme_light() +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Treatment") + ylab("Shannon") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))


z <- ggplot(data = sea33_means_dfSha, aes(x = season, y = emmean, fill = season)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  scale_fill_manual(values = my_colors2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 1)) + 
  xlab("Season") + ylab("Shannon") + 
  guides(fill=FALSE) +
  facet_grid(.~day) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  #theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(color="black", size=10,face= "bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) +
  theme(text=element_text(family="Times New Roman"))


ggarrange(a, b, c, d, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow=2, font.label = list(family = "Times New Roman"))
ggarrange(w, x, y, z, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow=2, font.label = list(family = "Times New Roman"))

          