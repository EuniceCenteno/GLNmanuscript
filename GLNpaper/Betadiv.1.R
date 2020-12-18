##----        R-Script for piglet ileal microbiota NMDS plots and ellipses
#                         Beta diversity statistical analysis
#                     Autor: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020

setwd('~/Desktop/eunice/Glutamine/No.D0/metadata/')
###### DATA #####
library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(fdrtool)
library(ggpubr)
library(reshape2)

##### function #####
## Using as reference "Tutorial on drawing an NMDS using ggplot2 by Umer Zeeshan Ijaz"
# http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ggplot2/NMDS.R

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

##----Create a plot to visualyze the beta diversity------ 

###nmds otu
#otu_table <-"No0.OTU.otusubsample.list.txt"
otu_table <-"No0.OTU.otusubsample.txt"
#metadata <- "No0.JayJohnson2.design.txt" --> no order in the metadata
metadata <- "No0.JayJohnson2.txt"

otu.subsample2 <- read.table(otu_table, header = TRUE) 
otu_subsample <- otu.subsample2
rownames(otu_subsample) <- otu_subsample$Group # stores the sample id info as the rownames of the dataframe rather

meta <- read.table(file = metadata, sep = '\t', header = TRUE)
meta <- meta[meta$Group %in% rownames(otu_subsample),]
meta <- meta[,-c(1,6)] ###remove columns that are no neccesary for this analysis 

### separating the METAdata by day
# DAY 13 
meta_d13 <- subset(meta, day == 13)
str(meta_d13)
meta_d13$season <- factor(meta_d13$season)
meta_d13$day <- factor(meta_d13$day)
meta_d13$Group <- factor(meta_d13$Group)
meta_d13$exp_season <- factor(meta_d13$exp_season)
meta_d13$exp_treatment <- factor(meta_d13$exp_treatment)
meta_d13$treatment <- factor(meta_d13$treatment)
levels(meta_d13$treatment) <- list("AB"="A","GLN"="GLN", "NA"="NA_C")
levels(meta_d13$treatment)
str(meta_d13)

#### prepating OTUs table d13 for NMDS plot
# this makes sure there are no samples in the OTU table that are not in our metadata
otu_subsampled13 <- otu_subsample[rownames(otu_subsample) %in% meta_d13$Group,]
str(otu_subsampled13)
otu_subsampled13$Group <- factor(otu_subsampled13$Group)
otu_subsampled13 <- otu_subsampled13[,-c(1:5)]  # removes extra info that mothur includes in their OTU tables

# DAY 33
meta_d33 <- subset(meta, day == 33)
str(meta_d33)
meta_d33$season <- factor(meta_d33$season)
meta_d33$day <- factor(meta_d33$day)
meta_d33$Group <- factor(meta_d33$Group)
meta_d33$exp_season <- factor(meta_d33$exp_season)
meta_d33$exp_treatment <- factor(meta_d33$exp_treatment)
meta_d33$treatment <- factor(meta_d33$treatment)
levels(meta_d33$treatment) <- list("AB"="A","GLN"="GLN", "NA"="NA_C")
str(meta_d33)

#### prepating OTUs table d33 for NMDS
# this makes sure there are no samples in the OTU table that are not in our metadata
otu_subsampled33 <- otu_subsample[rownames(otu_subsample) %in% meta_d33$Group,]
str(otu_subsampled33)
otu_subsampled33$Group <- factor(otu_subsampled33$Group)
otu_subsampled33 <- otu_subsampled33[,-c(1:5)]  # removes extra info that mothur includes in their OTU tables

##################################################
##################################################

# ----Calculation of the Distance Matrix ------
# Using Bray-Curtis distances for NMDS

# Day 13 Creating distance mitrix into a dataframe
dist.matr.brayd13 <- vegdist(otu_subsampled13, method = 'bray')
otu_subsample_metad13 <- merge(x = meta_d13, y = otu_subsampled13, by.x = "Group", by.y = 0)
str(otu_subsample_metad13)


# Day 33
dist.matr.brayd33 <- vegdist(otu_subsampled33, method = 'bray')
otu_subsample_metad33 <- merge(x = meta_d33, y = otu_subsampled33, by.x = "Group", by.y = 0)
str(otu_subsample_metad33)

#################
# Creating the NMDS ordination (Using the vegan function metaMDS)
# This helps to represent in two dimensions the distances or similatiries previously calculated using the distance matrix for each of the samples
mdsd13 <- metaMDS(dist.matr.brayd13, k = 2,trymax = 1000, autotransform = FALSE)
mdsd33 <- metaMDS(dist.matr.brayd33, k = 2,trymax = 1000, autotransform = FALSE)

# the stress of an ordination tells you how well represented are the distances or similarities from the distance matrix.
# You should expect to see a low stress value (<0.2) which indicate that the spatial distances in the ordination accuratelt represente the distances you calculated
mdsd13$stress
mdsd33$stress

#D13
nmdsd13 <-as.data.frame(mdsd13$points)
nmdsd13$Group <- rownames(nmdsd13) ########sample coordinates
metanmdsd13 <- merge(nmdsd13, otu_subsample_metad13, by.x = 'Group', by.y = 'Group')
metanmdsd13$day <- factor(metanmdsd13$day)
metanmdsd13$season <- factor(metanmdsd13$season)
metanmdsd13$treatment <- factor(metanmdsd13$treatment)
metanmdsd13$Group <- factor(metanmdsd13$Group)
str(metanmdsd13)

#D33
nmdsd33 <-as.data.frame(mdsd33$points)
nmdsd33$Group <- rownames(nmdsd33) ########sample coordinates
metanmdsd33 <- merge(nmdsd33, otu_subsample_metad33, by.x = 'Group', by.y = 'Group')
metanmdsd33$day <- factor(metanmdsd33$day)
metanmdsd33$season <- factor(metanmdsd33$season)
metanmdsd33$treatment <- factor(metanmdsd33$treatment)
metanmdsd33$Group <- factor(metanmdsd33$Group)
str(metanmdsd33)

# this generates a dataframe containing the group centroids for the ellipses
#####------CENTRIOIDS---------############

#D13 centroids for treatment and season
NMDS.meantd13 <- aggregate(metanmdsd13[,2:3], list(Group=metanmdsd13$treatment), mean)
colnames(NMDS.meantd13) <- c('treatmentC','grouptX', 'grouptY')

NMDS.meansd13 <- aggregate(metanmdsd13[,2:3], list(Group=metanmdsd13$season), mean)
colnames(NMDS.meansd13) <- c('seasonC','groupsX', 'groupsY')

# merging the group centroids with the rest of the NMDS data #
metanmdsd13 <- merge(NMDS.meantd13, metanmdsd13, by.x = "treatmentC", by.y= "treatment")
str(metanmdsd13)

metanmdsd13 <- merge(NMDS.meansd13, metanmdsd13, by.x = "seasonC", by.y= "season")
str(metanmdsd13)

#D33 centroids for treatment and season
NMDS.meantd33 <- aggregate(metanmdsd33[,2:3], list(Group=metanmdsd33$treatment), mean)
colnames(NMDS.meantd33) <- c('treatmentC','grouptX', 'grouptY')

NMDS.meansd33 <- aggregate(metanmdsd33[,2:3], list(Group=metanmdsd33$season), mean)
colnames(NMDS.meansd33) <- c('seasonC','groupsX', 'groupsY')

# merging the group centroids with the rest of the NMDS data #
metanmdsd33 <- merge(NMDS.meantd33, metanmdsd33, by.x = "treatmentC", by.y= "treatment")
str(metanmdsd33)

metanmdsd33 <- merge(NMDS.meansd33, metanmdsd33, by.x = "seasonC", by.y= "season")
str(metanmdsd33)

## Plotting each of the points using ggplot
ggplot(metanmdsd13, aes(x=MDS1, y=MDS2)) + geom_point(color='green')
ggplot(metanmdsd33, aes(x=MDS1, y=MDS2)) + geom_point(color='green')

ggplot(metanmdsd13, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=day))
ggplot(metanmdsd33, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=day))


########## ELLIPSE #############
# Check if the names in your metadata are in order
# this is how I check this:

nmdsd13$Group == metanmdsd13$Group  # by looking at these two metadata, we can check if the names are in order
# If the names are not in the sames order, we can run the next lines
metanmdsd13 <- metanmdsd13[match(nmdsd13$Group,metanmdsd13$Group),] 
nmdsd13$Group == metanmdsd13$Group  

nmdsd33$Group == metanmdsd33$Group  
metanmdsd33 <- metanmdsd33[match(nmdsd33$Group,metanmdsd33$Group),] 
nmdsd33$Group == metanmdsd33$Group 

# this generates a dataframe containing the group centroids
#  D13 
ordt13 <- ordiellipse(mdsd13, metanmdsd13$treatmentC, label = TRUE, conf = .95, kind = 'se', draw = 'none')
ords13 <- ordiellipse(mdsd13, metanmdsd13$seasonC, label = TRUE, conf = .95, kind = 'se', draw = 'none')


#  D33 
ordt33 <- ordiellipse(mdsd33, metanmdsd33$treatmentC, label = TRUE, conf = .95, kind = 'se', draw = 'none')
ords33 <- ordiellipse(mdsd33, metanmdsd33$seasonC, label = TRUE, conf = .95, kind = 'se', draw = 'none')

#########-------ellipse coordinates--------######

#Day 13
#- treatment
df_d13T <- data.frame()
for (d in levels(metanmdsd13$treatmentC)){
  df_d13T <- rbind(df_d13T, cbind(as.data.frame(with(metanmdsd13[metanmdsd13$treatmentC == d,],
                                                     veganCovEllipse(ordt13[[d]]$cov, ordt13[[d]]$center, ordt13[[d]]$scale))), treatmentC=d))
}

dTR_d13 <- merge(df_d13T, metanmdsd13, by.x = 'treatmentC', by.y='treatmentC')
dTR_d13 <- dTR_d13[,-c(4:6, 14:1021)] ###remove the OTUs
colnames(dTR_d13) <- c('treatment', 'tEMDS1', 'tEMDS2','groupTX', 'groupTY', 'Group', 'NMDS1', 'NMDS2', 'day', 'exp_group' ) # just making it so our column names are consistent

# season
df_d13S <- data.frame()
for (d in levels(metanmdsd13$seasonC)){
  df_d13S <- rbind(df_d13S, cbind(as.data.frame(with(metanmdsd13[metanmdsd13$seasonC == d,],
                                                     veganCovEllipse(ords13[[d]]$cov, ords13[[d]]$center, ords13[[d]]$scale))), seasonC=d))
}

dS_d13 <- merge(df_d13S, metanmdsd13, by.x = 'seasonC', by.y='seasonC')
dS_d13 <- dS_d13[,-c(6:8,14:1021)] ###remove the OTUs
colnames(dS_d13) <- c('season', 'sEMDS1', 'sEMDS2', 'groupSX', 'groupSY', 'Group', 'NMDS1', 'NMDS2', 'day', 'exp_group' ) # just making it so our column names are consistent

##Day 33
df_d33T <- data.frame()
for (d in levels(metanmdsd33$treatmentC)){
  df_d33T <- rbind(df_d33T, cbind(as.data.frame(with(metanmdsd33[metanmdsd33$treatmentC == d,],
                                                     veganCovEllipse(ordt33[[d]]$cov, ordt33[[d]]$center, ordt33[[d]]$scale))), treatmentC=d))
}

dTR_d33 <- merge(df_d33T, metanmdsd33, by.x = 'treatmentC', by.y='treatmentC')
dTR_d33 <- dTR_d33[,-c(4:6, 14:1021)] ###remove the OTUs
colnames(dTR_d33) <- c('treatment', 'tEMDS1', 'tEMDS2','groupTX', 'groupTY', 'Group', 'NMDS1', 'NMDS2', 'day', 'exp_group' ) # just making it so our column names are consistent

# season
df_d33S <- data.frame()
for (d in levels(metanmdsd33$seasonC)){
  df_d33S <- rbind(df_d33S, cbind(as.data.frame(with(metanmdsd33[metanmdsd33$seasonC == d,],
                                                     veganCovEllipse(ords33[[d]]$cov, ords33[[d]]$center, ords33[[d]]$scale))), seasonC=d))
}

dS_d33 <- merge(df_d33S, metanmdsd33, by.x = 'seasonC', by.y='seasonC')
dS_d33 <- dS_d33[,-c(6:8,14:1021)] ###remove the OTUs
colnames(dS_d33) <- c('season', 'sEMDS1', 'sEMDS2','groupSX', 'groupSY', 'Group', 'NMDS1', 'NMDS2', 'day', 'exp_group' ) # just making it so our column names are consistent

str(dTR_d13)
str(dS_d13)
str(dTR_d33)
str(dS_d33)

#####EMDS1 and EMDS2 are the ellipse coordinates
#####groupX and groupY are the ellipse centroids
#####NMDS1 and NMDS2 are the coordinates for the samples

################ellipse using geom_polygon
#####EMDS1 and EMDS2 are the ellipse coordinates
#####groupX and groupY are the ellipse centroids
#####NMDS1 and NMDS2 are the coordinates for the samples

###################ellipse using stat_ellipse
##------PLOTTING ELLIPSES
str(dTR_d13)
dTR_d13$treatment <- factor(dTR_d13$treatment)
str(dS_d13$season)
dS_d13$season <- factor(dS_d13$season)
levels(dS_d13$season) <- list("Summer"="Summer", "Spring"="Spring")
str(dTR_d33)
dTR_d33$treatment <- factor(dTR_d33$treatment)
str(dS_d33$season)
dS_d33$season <- factor(dS_d33$season)
levels(dS_d33$season) <- list("Summer"="Summer", "Spring"="Spring")

## DAY 13, TREATMENT
library(extrafont)
font_import()
loadfonts(device = "win")
fonts()
my_colors <- c('yellowgreen',"sienna2", 'skyblue1') 

a <-ggplot(metanmdsd13, aes(x=MDS1, y=MDS2)) +
  theme_bw()+
  geom_point(data = dTR_d13, aes(x=NMDS1, y= NMDS2,color=treatment), size =1) +
  scale_color_manual(values = my_colors) +
  geom_segment(data= dTR_d13, aes(x=NMDS1, y=NMDS2, xend=groupTX, yend=groupTY, color= treatment), size = .05) + 
  stat_ellipse(geom = "polygon", data = dTR_d13, aes(x=tEMDS1, y=tEMDS2, group=treatment, fill= treatment), alpha = 0.2) + 
  scale_fill_manual(values = my_colors) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mdsd13$stress, digits = 2))) +
  labs(color= "Treatment") +
  facet_grid(.~day) +
  guides(fill=FALSE) +
  #guides(color=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  theme(text=element_text(family="Times New Roman"))


#DAY 13 season
my_colors2 <- c('tan3',"turquoise4", 'skyblue1') 
c <-ggplot(metanmdsd13, aes(x=MDS1, y=MDS2)) +
  theme_bw()+
  geom_point(data = dS_d13, aes(x=NMDS1, y= NMDS2,color=season), size =1) +
  scale_color_manual(values = my_colors2) +
  geom_segment(data= dS_d13, aes(x=NMDS1, y=NMDS2, xend=groupSX, yend=groupSY, color= season), size = .05) + 
  stat_ellipse(geom = "polygon", data = dS_d13, aes(x=sEMDS1, y=sEMDS2, group=season, fill= season), alpha = 0.2) + 
  scale_fill_manual(values = my_colors2) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mdsd13$stress, digits = 2))) +
  labs(color= "Season") +
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) +
  theme(text=element_text(family="Times New Roman"))

## DAY 33, TREATMENT
b <- ggplot(metanmdsd33, aes(x=MDS1, y=MDS2)) +
  theme_bw()+
  geom_point(data = dTR_d33, aes(x=NMDS1, y= NMDS2,color=treatment), size =1) +
  scale_color_manual(values = my_colors) +
  geom_segment(data= dTR_d33, aes(x=NMDS1, y=NMDS2, xend=groupTX, yend=groupTY, color= treatment), size = .05) + 
  stat_ellipse(geom = "polygon", data = dTR_d33, aes(x=tEMDS1, y=tEMDS2, group=treatment, fill= treatment), alpha = 0.2) + 
  scale_fill_manual(values = my_colors) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mdsd33$stress, digits = 2))) +
  labs(color= "Treatment") +
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) + 
  theme(text=element_text(family="Times New Roman"))

#DAY 33 season
d <- ggplot(metanmdsd33, aes(x=MDS1, y=MDS2)) +
  theme_bw()+
  geom_point(data = dS_d33, aes(x=NMDS1, y= NMDS2,color=season), size =1) +
  scale_color_manual(values = my_colors2) +
  geom_segment(data= dS_d33, aes(x=NMDS1, y=NMDS2, xend=groupSX, yend=groupSY, color= season), size = .05) + 
  stat_ellipse(geom = "polygon", data = dS_d33, aes(x=sEMDS1, y=sEMDS2, group=season, fill= season), alpha = 0.2) + 
  scale_fill_manual(values = my_colors2) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mdsd33$stress, digits = 2))) +
  labs(color= "Season") +
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 10)) +
  theme(text=element_text(family="Times New Roman"))


ggarrange(a, b, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))
ggarrange(c, d, labels = c("A", "B"),
          ncol = 2, font.label = list(family = "Times New Roman"))

####### Distance calculates with pairwsie comparison
dist.t13 <- "dist.13.2.txt"
dist.T13 <- read.table(file = dist.t13, sep = '\t', header = TRUE)
str(dist.T13)
dist.T13$Group1 <- factor(dist.T13$Group1)
dist.T13$Group2 <- factor(dist.T13$Group2)
dist.T13$group <- factor(dist.T13$group)
dist.T13$pairwise <- factor(dist.T13$pairwise)
levels(dist.T13$pairwise)
dist.T13 <- dist.T13[,-c(1:2)]
levels(dist.T13$pairwise) <- list("AB-AB"="AB-AB", "GLN-GLN"="GLN-GLN", "NA-NA"="NA-NA", "AB-NA"="AB-NA", "GLN-AB"="GLN-AB", "GLN-NA"="GLN-NA")
dist.T13$day <- 13
dist.T13$day <- factor(dist.T13$day)
str(dist.T13)
my_colors3 <- c("gray0", "gray24", "gray43", 'gray65', "gray86", "gray95") 

SuT <- ggplot(dist.T13, aes(pairwise, value, fill=pairwise))+
  stat_boxplot( aes(pairwise, value), 
                geom='errorbar', linetype=1, width=0.5) + #whiskers
  geom_boxplot( aes(pairwise, value),outlier.shape=1) +
  theme_light() +
  scale_fill_manual(values= my_colors3) +
  xlab("Treatment") + ylab("Distance") + 
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  theme(text=element_text(family="Times New Roman"))

dist.t13 <- "dist.t33.txt"
dist.T33 <- read.table(file = dist.t13, sep = '\t', header = TRUE)
str(dist.T33)
dist.T33$Group1 <- factor(dist.T33$Group1)
dist.T33$Group2 <- factor(dist.T33$Group2)
dist.T33$pairwise <- factor(dist.T33$pairwise)
levels(dist.T33$pairwise)
dist.T33 <- dist.T33[,-c(1:2)]
levels(dist.T33$pairwise) <- list("AB-AB"="AB-AB", "GLN-GLN"="GLN-GLN", "NA-NA"="NA-NA", "AB-NA"="AB-NA", "GLN-AB"="GLN-AB", "GLN-NA"="GLN-NA")
dist.T33$day <- 33
dist.T33$day <- factor(dist.T33$day)
str(dist.T33)
levels(dist.T33$pairwise)
is.na(dist.T33)

 
SpT <- ggplot(dist.T33, aes(pairwise, value, fill=pairwise))+
  stat_boxplot( aes(pairwise, value), 
                geom='errorbar', linetype=1, width=0.5) + #whiskers
  geom_boxplot( aes(pairwise, value),outlier.shape=1) +
  theme_light() +
  scale_fill_manual(values = my_colors3) +
  xlab("Treatment") + ylab("Distance") + 
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 11, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 11, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 11, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  theme(text=element_text(family="Times New Roman"))

ggarrange(SuT, SpT,labels = c("A", "B"),
          nrow = 2, font.label = list(family = "Times New Roman"))

## Season
dist.s13 <- "dist.s13.2.txt"
dist.S13 <- read.table(file = dist.s13, sep = '\t', header = TRUE)
str(dist.S13)
dist.S13$Group1 <- factor(dist.S13$Group1)
dist.S13$Group2 <- factor(dist.S13$Group2)
dist.S13$pairwise <- factor(dist.S13$pairwise)
levels(dist.S13$pairwise)
dist.S13 <- dist.S13[,-c(1:2)]
levels(dist.S13$pairwise) <- list("Sum-Sum"="Sum-Sum", "Spr-Spr"="Spr-Spr", "Sum-Spr"="Sum-Spr")
dist.S13$day <- 13
dist.S13$day <- factor(dist.S13$day)
str(dist.S13)

Su1 <- ggplot(dist.S13, aes(pairwise, value, fill=pairwise))+
  stat_boxplot( aes(pairwise, value), 
                geom='errorbar', linetype=1, width=0.5) + #whiskers
  geom_boxplot( aes(pairwise, value),outlier.shape=1) +
  theme_light() +
  scale_fill_manual(values = c("black", "gray", "white")) +
  xlab("Treatment") + ylab("Distance") + 
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  theme(text=element_text(family="Times New Roman"))

dist.s33 <- "dist.s33.2.txt"
dist.S33 <- read.table(file = dist.s33, sep = '\t', header = TRUE)
str(dist.S33)
dist.S33$Group1 <- factor(dist.S33$Group1)
dist.S33$Group2 <- factor(dist.S33$Group2)
dist.S33$pairwise <- factor(dist.S33$pairwise)
levels(dist.S33$pairwise)
dist.S33 <- dist.S33[,-c(1:2)]
levels(dist.S33$pairwise) <- list("Sum-Sum"="Sum-Sum", "Spr-Spr"="Spri-Spr", "Sum-Spr"="Sum-Spr")
dist.S33$day <- 33
dist.S33$day <- factor(dist.S33$day)
str(dist.S33)

Su2 <- ggplot(dist.S33, aes(pairwise, value, fill=pairwise))+
  stat_boxplot( aes(pairwise, value), 
                geom='errorbar', linetype=1, width=0.5) + #whiskers
  geom_boxplot( aes(pairwise, value),outlier.shape=1) +
  theme_light() +
  scale_fill_manual(values = c("black", "gray","white")) +
  xlab("Treatment") + ylab("Distance") + 
  facet_grid(.~day) +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10, face="bold"), axis.text.y = element_text(color = "black", size = 10)) +
  theme(axis.title.x = element_blank(), axis.ticks = element_blank()) + 
  theme(text=element_text(family="Times New Roman"))

ggarrange(Su1, Su2,labels = c("A", "B"),
          nrow = 2, font.label = list(family = "Times New Roman"))

#####-----STATISTICAL ANALYSIS---------
## PERMANOVA using adonis function

#subet beta diversity for day 13 
meta %>% tally()
meta %>% count(day)
# 13       56
# 33       51

## We need the metadata and the OTU table
str(meta_d13)
levels(meta_d13$treatment) <- list("A"="A", "GLN"="GLN", "NA"="NA")
levels(meta_d13$season) <- list("Summer"="Summer", "Spring"="Spring")
str(otu_subsampled13)

#PERMANOVA 
Per_d13 <- adonis(otu_subsampled13 ~ meta_d13$treatment + meta_d13$season + meta_d13$treatment:meta_d13$season, permutations = 999)
Per_d13

#we can also performe a pairwise comparison with the function 
# Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
#https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

# Pairwise comparison
comp.trt.D13 <- pairwise.adonis(otu_subsampled13, meta_d13$treatment, sim.method = "bray", p.adjust.m = "BH")
comp.trt.D13

comp.season.D13 <- pairwise.adonis(otu_subsampled13, meta_d13$season, sim.method = "bray", p.adjust.m = "BH")
comp.season.D13

comp.int.D13 <- pairwise.adonis(otu_subsampled13,  meta_d13$treatment:meta_d13$season, sim.method = "bray", p.adjust.m = "BH")
comp.int.D13

#subset beta diversity for day 33
# we need the metadata and the OTU table
str(meta_d13)
str(meta_d33)
levels(meta_d33$treatment) <- list("A"="A", "GLN"="GLN", "NA"="NA")
levels(meta_d33$season) <- list("Summer"="Summer", "Spring"="Spring")
levels(meta_d33$treatment)
str(otu_subsampled33)

#PERMANOVA 
Per_d33 <- adonis(otu_subsampled33 ~ meta_d33$treatment + meta_d33$season + meta_d33$treatment:meta_d33$season, permutations = 999)
Per_d33

### Pairwise comparison
comp.trt.D33 <- pairwise.adonis(otu_subsampled33, meta_d33$treatment, sim.method = "bray", p.adjust.m = "BH")
comp.trt.D33

comp.season.D33 <- pairwise.adonis(otu_subsampled33, meta_d33$season, sim.method = "bray", p.adjust.m = "BH")
comp.season.D33

####---------------------------- DISPERSION 
### Using vegan's betadisper() function ###

#check homogeneity of variance
#------D13- TREATMENT
dispers13T <- betadisper(dist.matr.brayd13, type = c("centroid"), group = meta_d13$treatment)
dispers13T
boxplot(dispers13T)

dispers13T.2 <- data.frame(Distance_to_centroid=dispers13T$distances,Group=dispers13T$group)
groups <-dispers13T$group

pdispers13T <- permutest(dispers13T, permutations = 999)#, pairwise = TRUE)
pdispers13T

#------D13- SEASON
levels(meta_d13$season)
levels(meta_d13$season) <- list("Summer"="Summer", "Spring"="Spring")
dispers13S <- betadisper(dist.matr.brayd13,type = c("centroid"), group = meta_d13$season)
dispers13S
boxplot(dispers13S)

pdispers13S <- permutest(dispers13S, permutations = 999)#, pairwise = TRUE)
pdispers13S

#------D33- TREATMENT
dispers33T <- betadisper(dist.matr.brayd33,type = c("centroid"),  group = meta_d33$treatment)
dispers33T
boxplot(dispers33T)

pdispers33T <- permutest(dispers33T, permutations = 999)#, pairwise = TRUE)
pdispers33T

#------D33- SEASON
levels(meta_d33$season)
levels(meta_d33$season) <- list("Summer"="Summer", "Spring"="Spring")
dispers33S <- betadisper(dist.matr.brayd33, type = c("centroid"), group = meta_d33$season)
dispers33S
boxplot(dispers33S)

pdispers33S <- permutest(dispers33S, permutations = 999)#, pairwise = TRUE)
pdispers33S
