###               R-Script for Piglet ileal microbiota Taxa Barplot
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020
# rcenteno@purdue.edu
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)

setwd("~/Desktop/eu/TesisGLN/new.glutamine/glutamine.list/No.D0/metadata/")
##

##### Taxonomy barplot
Otufile <- "No0.PHY.otusubsample.list.txt"
taxfile <- "PHY.tax.list.txt"
metadata <- "No0.JayJohnson2.design.txt"

#first the data file
otu_subsample <- read.table(Otufile, sep = "\t", header = T) 
otu_subsample <- otu_subsample[,c(-1,-3)]
row.names(otu_subsample) <- otu_subsample$Group
otu_subsample <- otu_subsample[,-1]
#there are 246 OTUs
#then the taxonomy file
taxonomy <- read.table(taxfile, sep = "\t", header = T) 
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

#then the metadata file
meta <- read.table(file = metadata, sep = '\t', header = TRUE)
meta <- arrange(meta, day, treatment, season)
meta$treatment <- factor(meta$treatment, levels = c("A", "GLN", "NA_C"))
str(meta)
levels(meta$treatment) <- list("A"="A", "GLN"="GLN", "NA"="NA_C")
order_groups <- meta$group
meta$groups <- factor(meta$group, levels = order_groups)

### CALCULATION OF THE ABUNDANCE OF EACH OTU  
otu.summary <- prop.table(as.matrix(otu_subsample), 1) 
otu_abund <- colSums(otu.summary)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
melt_otu <- melt(otu.summary_sorted[,c(1:246)]) ###TOTAL NUMBER OF OTUS
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")

#merging the abundance of each OTU with the metadata and the taxonomy file
meta_otu <- merge(meta, melt_otu, by.x = "Group", by.y = "Sample")
meta_otu <- meta_otu[-c(8)]
meta_otu_tax <- merge(meta_otu, taxonomy)
str(meta_otu_tax)
meta_otu_tax<- meta_otu_tax[-c(17)]
meta_otu_tax$treatment <- factor(meta_otu_tax$treatment, levels = c("A", "GLN", "NA"))
meta_otu_tax$group <- factor(meta_otu_tax$group, levels = order_groups)
summary(meta_otu_tax$group) ###to check that all the samples have the same number of OTUs (246 total) 
meta_otu_tax$season <- factor(meta_otu_tax$season)
str(meta_otu_tax)
levels(meta_otu_tax$season) <- list("Summer"="Summer", "Spring"="Spring")
str(meta_otu_tax)
meta_otu_tax$family <- factor(meta_otu_tax$family)
meta_otu_tax$genus <- factor(meta_otu_tax$genus)

###----- TAXONOMY CALCULATION FOR DAY 13

##this is to calculate the abundance of each of the phyla and family groups on day 13
meta_otu_13 <- subset(meta_otu_tax, day== "13")
str(meta_otu_13)
meta_otu_13$group <- factor(meta_otu_13$group)
meta_otu_13$OTU <- factor(meta_otu_13$OTU)
meta_otu_13$family <- factor(meta_otu_13$family)
meta_otu_13$genus <- factor(meta_otu_13$genus)
meta_otu_13$season <- factor(meta_otu_13$season)
meta_otu_13$treatment <- factor(meta_otu_13$treatment)
meta_otu_13$day <- factor(meta_otu_13$day)

## The abundance at a family level
str(meta_otu_13)
family13 <- meta_otu_13 %>% 
  group_by(Group, family) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(family) %>%
  summarise(taxa.average = mean(taxa.sum)) 
attach(family13)
family13 <- family13[order(-taxa.average),]
#write.table(family13, file = "familylevel.txt", sep = "\t", quote = FALSE)

### Abundance at a phlyum level
phylum13 <- meta_otu_13 %>% 
  group_by(Group, phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(phylum) %>%
  summarise(taxa.average = mean(taxa.sum)) 
attach(phylum13)
phylum13 <- phylum13[order(-taxa.average),]
#write.table(phylum13, file = "phylumlevel13.txt", sep = "\t", quote = FALSE)

### Abundance at a genus level
genus13 <- meta_otu_13 %>% 
  group_by(genus) %>% 
  summarise(Abundance.average = mean(Abundance))
attach(genus13)
genus13 <- genus13[order(-Abundance.average),]
#write.table(genus13, file = "genuslevel13.txt", sep = "\t", quote = FALSE)

# PLOT FOR THE FIRST 10 GENUS
genus13.top <- genus13[1:10,]
str(genus13.top)
genus13.top$genus <- factor(genus13.top$genus)
levels(genus13.top$genus)
levels(genus13.top$genus) <- list("Acidaminococcus"="Acidaminococcus(100)", "Actinobacillus"="Actinobacillus(100)",
                                     "Bifidobacterium"="Bifidobacterium(100)","Clostridium sensu stricto"="Clostridium_sensu_stricto(100)", "Enterobacteriaceae unclass"="Enterobacteriaceae_unclassified(100)",
                                     "Lactobacillus"="Lactobacillus(100)","Megasphaera"="Megasphaera(100)", "Mitsuokella"="Mitsuokella(100)", "Olsenella"="Olsenella(100)",
                                     "Pasteurellaceae unclass"="Pasteurellaceae_unclassified(100)", "Prevotella"="Prevotella(100)", "Streptococcus"="Streptococcus(100)","Terrisporobacter"="Terrisporobacter(100)", 
                                     "Veillonella"="Veillonella(100)", "Veillonellaceae unclass"="Veillonellaceae_unclassified(100)")

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

## Treatment effect (Top 10) day 13
ggplot(genus13.top, aes(x = 0, y = Abundance.average, fill = genus)) + 
  theme_bw()+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(fill= "Genus") +
  ylab("Relative Abundance") + xlab("Ileal Bacterial Composition") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10,face= "bold")) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=10,face= "bold")) 


### Barplot for the top 10 "family group"
str(family13)
family13.top <- family13[1:10,]
family13.top$family <- factor(family13.top$family)
levels(family13.top$family)
levels(family13.top$family) <- list("Bifidobacteriaceae"="Bifidobacteriaceae(100)", "Clostridiaceae 1"="Clostridiaceae_1(100)", "Coriobacteriaceae"="Coriobacteriaceae(100)",
                                     "Enterobacteriaceae"="Enterobacteriaceae(100)","Lactobacillaceae"="Lactobacillaceae(100)",
                                     "Pasteurellaceae"="Pasteurellaceae(100)", "Peptostreptococcaceae"="Peptostreptococcaceae(100)", "Prevotellaceae"="Prevotellaceae(100)",
                                     "Streptococcaceae"="Streptococcaceae(100)", "Veillonellaceae"="Veillonellaceae(100)")
my_colors2 <- c("palegreen3", "darkgoldenrod2", "thistle","#AD6F3B", "powderblue","lemonchiffon3", "slategray1","lightgoldenrodyellow",
                "darkseagreen3", "ivory3", "lightpink4", "mistyrose3", "salmon3", "pink")    

ggplot(family13.top, aes(x = 1, y = taxa.average, fill = family)) + 
  theme_bw()+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors2) +
  labs(fill= "Family") +
  ylab("Relative Abundance") + xlab("Ileal Bacterial Composition") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10,face= "bold")) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=10,face= "bold")) 


###----- TAXONOMY CALCULATION FOR DAY 33

##this is to calculate the abundance of each of the phyla and family groups on day 33
meta_otu_33 <- subset(meta_otu_tax, day== "33")
str(meta_otu_33)
meta_otu_33$group <- factor(meta_otu_33$group)
meta_otu_33$OTU <- factor(meta_otu_33$OTU)
meta_otu_33$family <- factor(meta_otu_33$family)
meta_otu_33$genus <- factor(meta_otu_33$genus)
meta_otu_33$season <- factor(meta_otu_33$season)
meta_otu_33$treatment <- factor(meta_otu_33$treatment)
meta_otu_33$day <- factor(meta_otu_33$day)

## The abundance at a family level
str(meta_otu_33)
family33 <- meta_otu_33 %>% 
  group_by(Group, family) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(family) %>%
  summarise(taxa.average = mean(taxa.sum)) 
attach(family33)
family33 <- family33[order(-taxa.average),]
#write.table(family13, file = "familylevel.txt", sep = "\t", quote = FALSE)

### Abundance at a phlyum level
phylum33 <- meta_otu_33 %>% 
  group_by(Group, phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(phylum) %>%
  summarise(taxa.average = mean(taxa.sum)) 
attach(phylum33)
phylum33 <- phylum33[order(-taxa.average),]
#write.table(phylum13, file = "phylumlevel13.txt", sep = "\t", quote = FALSE)

### Abundance at a genus level
genus33 <- meta_otu_33 %>% 
  group_by(genus) %>% 
  summarise(Abundance.average = mean(Abundance))
attach(genus33)
genus33 <- genus33[order(-Abundance.average),]
#write.table(genus13, file = "genuslevel13.txt", sep = "\t", quote = FALSE)

# PLOT FOR THE FIRST 10 GENUS
genus33.top <- genus33[1:10,]
str(genus33.top)
genus33.top$genus <- factor(genus33.top$genus)
levels(genus33.top$genus)
levels(genus33.top$genus) <- list("Actinobacillus"="Actinobacillus(100)","Clostridium sensu stricto"="Clostridium_sensu_stricto(100)", "Enterobacteriaceae unclass"="Enterobacteriaceae_unclassified(100)",
                                  "Lactobacillus"="Lactobacillus(100)","Megasphaera"="Megasphaera(100)", "Pasteurellaceae unclass"="Pasteurellaceae_unclassified(100)", "Prevotella"="Prevotella(100)", 
                                  "Streptococcus"="Streptococcus(100)","Terrisporobacter"="Terrisporobacter(100)", "Veillonella"="Veillonella(100)")

## Treatment effect (Top 10) day 13
ggplot(genus33.top, aes(x = 1, y = Abundance.average, fill = genus)) + 
  theme_bw()+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  labs(fill= "Genus") +
  ylab("Relative Abundance") + xlab("Ileal Bacterial Composition") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10,face= "bold")) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=10,face= "bold")) 


### Barplot for the top 10 "family group"
str(family33)
family33.top <- family33[1:10,]
family33.top$family <- factor(family33.top$family)
levels(family33.top$family)
levels(family33.top$family) <- list("Bifidobacteriaceae"="Bifidobacteriaceae(100)", "Clostridiaceae 1"="Clostridiaceae_1(100)",
                                    "Enterobacteriaceae"="Enterobacteriaceae(100)", "Helicobacteraceae"="Helicobacteraceae(100)", "Lactobacillaceae"="Lactobacillaceae(100)",
                                    "Pasteurellaceae"="Pasteurellaceae(100)", "Peptostreptococcaceae"="Peptostreptococcaceae(100)", "Prevotellaceae"="Prevotellaceae(100)",
                                    "Streptococcaceae"="Streptococcaceae(100)", "Veillonellaceae"="Veillonellaceae(100)")

## Top 10 family
ggplot(family33.top, aes(x = 1, y = taxa.average, fill = family)) + 
  theme_bw()+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors2) +
  labs(fill= "Family") +
  ylab("Relative Abundance") + xlab("Ileal Bacterial Composition") + 
  theme(strip.text = element_text(size = 9, face = "bold", color= "black")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size = 10, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.y = element_text(color="black", size=10,face= "bold")) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_text(color="black", size=9,face= "bold")) 

