##----          R-Script for plotting the piglet ileal METASTATS results 
#                         Alpha diversity statistical analysis
#                     Author: Ruth Eunice Centeno Martinez, 2020
#          Written at Dr. Tim Johnson Lab, Dept. of Animal Sciences,       #
#                             Purdue University, 2020


library(ggplot2)

setwd('~/Desktop/eu/TesisGLN/new.glutamine/glutamine.list/No.D0/metasts.new/')

#### ---------- SEASON EFFECT----------------------
mets.season <- "metas.seaeffect.new.txt"
mets.season1 <- read.table(mets.season, sep = "\t", header = T) 
str(mets.season1)
mets.season1$day <- factor(mets.season1$day)

mets.season1$Genus <- factor(mets.season1$Genus)
mets.season1$season <- factor(mets.season1$season)
levels(mets.season1$season) <- list("Summer"="Summer", "Spring"="Spring")

str(mets.season1)

data <- data.frame(
  name=mets.season1$Genus,
  ave.abundance=mets.season1$adj.mean,
  sd=mets.season1$std,
  season=mets.season1$season,
  day=mets.season1$day
)

##subset for day 13
data_s13 <- subset(data, day == 13 )
str(data_s13)
data_s13$name <- factor(data_s13$name)
data_s13$day <- factor(data_s13$day)

### command that I used to order the genera from the one with high abundance to low abundace
levels(data_s13$name) <- list("Megasphaera"="Megasphaera", "Terrisporobacter"="Terrisporobacter","Actinobacillus"="Actinobacillus", "Veillonellaceae unclass"="Veillonellaceae unclass", "Acidaminococcus"="Acidaminococcus",
                              "Lachnospiraceae unclass"="Lachnospiraceae unclass", "Dialister"="Dialister", "Sharpea"="Sharpea", "Neisseria"="Neisseria", "Fusobacterium"="Fusobacterium", "Moraxella"="Moraxella",
                              "Aeromonas"="Aeromonas", "Treponema"="Treponema", "Pseudoramibacter"="Pseudoramibacter","Alloiococcus"="Alloiococcus", "Lactococcus"="Lactococcus", "Coriobacteriaceae unclass"="Coriobacteriaceae unclass",
                              "Allisonella"="Allisonella", "Succiniclasticum"= "Succiniclasticum", "Bacteroidetes unclass"="Bacteroidetes unclass", "Fusicatenibacter"="Fusicatenibacter","Weissella"="Weissella", "Anaerostipes"="Anaerostipes", "Intestinimonas"="Intestinimonas",
                              "Turicibacter"="Turicibacter", "Subdoligranulum"="Subdoligranulum" )

data_s13$season <- as.factor(data_s13$season)
data_s13$day <- as.factor(data_s13$day)

##season effect d13
# ggplot for all the species
ggplot(data = data_s13, aes(x=name, y=ave.abundance,  fill=season)) +
  geom_bar(aes(x=name, y=ave.abundance), position= position_dodge(1), stat="identity",  alpha=0.7) +
  theme_test() +
  scale_color_manual(values=c("#8b58c8","#aeb646"))+
  scale_fill_manual(values=c("#8b58c8","#aeb646"))+
  geom_errorbar( aes(x=name, ymin=ave.abundance-sd, ymax=ave.abundance+sd), position = position_dodge(1), stat = "identity", size=0.4, width=0.5) +
  theme(axis.text.x=element_text(colour="black")) + 
  scale_y_continuous() +
  theme(axis.text.y=element_text(colour="black")) + 
  xlab("Genus") + ylab("Relative Abundance") + 
  labs(fill= "Season") +
  facet_wrap(vars(name), scales = "free", strip.position = "bottom") +
  xlab(NULL) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.title = element_text(size = 10, face= "bold"))+
  theme(axis.text.x = element_text(color="black", size=8, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"))


###
##subset for day 33
data_s33 <- subset(data, day == 33 )
str(data_s33)
data_s33$name <- factor(data_s33$name)
data_s33$day <- factor(data_s33$day)

levels(data_s33$name) <- list("Lactobacillus"="Lactobacillus", "Terrisporobacter"="Terrisporobacter", "Veillonella"="Veillonella", "Olsenella"="Olsenella", "Mitsuokella"="Mitsuokella",
                              "Peptostreptococcaceae unclass"="Peptostreptococcaceae unclass", "Escherichia/Shigella"="Escherichia/Shigella", "Bifidobacterium"="Bifidobacterium","Ruminococcus"="Ruminococcus",
                              "Pseudoscardovia"="Pseudoscardovia","Planctomycetaceae unclass"="Planctomycetaceae unclass", "Collinsella"="Collinsella", "Stenotrophomonas"="Stenotrophomonas", "Streptobacillus"="Streptobacillus",
                              "Kitasatospora"="Kitasatospora", "Bradyrhizobium"="Bradyrhizobium", "Bacteroidetes unclass"="Bacteroidetes unclass", "Neisseria"="Neisseria","Bacteroides"="Bacteroides", "Selenomonas"="Selenomonas",
                              "Elusimicrobium"="Elusimicrobium","Brevundimonas"="Brevundimonas" )

data_s33$season <- as.factor(data_s33$season)
data_s33$day <- as.factor(data_s33$day)

##season effect d33
# ggplot for all the species
ggplot(data = data_s33, aes(x=name, y=ave.abundance,  fill=season)) +
  geom_bar(aes(x=name, y=ave.abundance), position= position_dodge(1), stat="identity",  alpha=0.7) +
  theme_test() +
  scale_color_manual(values=c("#8b58c8","#aeb646"))+
  scale_fill_manual(values=c("#8b58c8","#aeb646"))+
  geom_errorbar( aes(x=name, ymin=ave.abundance-sd, ymax=ave.abundance+sd), position = position_dodge(1), stat = "identity", size=0.4, width=0.5) +
  theme(axis.text.x=element_text(colour="black")) + 
  scale_y_continuous() +
  theme(axis.text.y=element_text(colour="black")) + 
  xlab("Genus") + ylab("Relative Abundance") + 
  labs(fill= "Season") +
  facet_wrap(vars(name), scales = "free", strip.position = "bottom") +
  xlab(NULL) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.title = element_text(size = 10, face= "bold"))+
  theme(axis.text.x = element_text(color="black", size=8, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"))

###------------------ TREATMENT EFFECT --------------

##treatment effect
mets.trt <- "metas.trteffect.new.txt"
mets.trt1 <- read.table(mets.trt, sep = "\t", header = T) 
str(mets.trt1)
levels(mets.trt1$treatment) <- list("A"="A", "GLN"="GLN", "NA"="NAS")
str(mets.trt1)
str(mets.trt1$Genus)
list(mets.trt1$Genus)

data2 <- data.frame(
  name=mets.trt1$Genus,
  ave.abundance=mets.trt1$adj.mean,
  sd=mets.trt1$std,
  treatment=mets.trt1$treatment,
  day=mets.trt1$day
)

##subset for 13
data_trt13 <- subset(data2, day == "13" )
str(data_trt13)
data_trt13$name <- as.factor(data_trt13$name)

levels(data_trt13$name) <- list("Lactobacillus"="Lactobacillus", "Helicobacter"="Helicobacter", "Leptotrichia"="Leptotrichia",
                                "Kitasatospora"="Kitasatospora","Actinomyces"="Actinomyces","Alloiococcus"="Alloiococcus","Treponema"="Treponema",
                                "Moraxella"="Moraxella","Clostridiales Incertae Sedis XIII unclass"="Clostridiales Incertae Sedis XIII unclass", "Fusicatenibacter"="Fusicatenibacter",
                                "Neisseriaceae unclass"="Neisseriaceae unclass")
str(data_trt13)
data_trt13$day <- as.factor(data_trt13$day)

my_colors <- c('yellowgreen',"sienna2", 'skyblue1') 

## ggplot for all the species
ggplot(data = data_trt13, aes(x=name, y=ave.abundance,  fill=treatment)) +
  geom_bar(aes(x=name, y=ave.abundance), position= position_dodge(1), stat="identity",  alpha=0.7) +
  theme_test() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  geom_errorbar( aes(x=name, ymin=ave.abundance-sd, ymax=ave.abundance+sd), position = position_dodge(1), stat = "identity", size=0.4, width=0.5) +
  theme(axis.text.x=element_text(colour="black")) + 
  scale_y_continuous() +
  theme(axis.text.y=element_text(colour="black")) + 
  xlab("Genus") + ylab("Relative Abundance") + 
  labs(fill= "Treatment") +
  facet_wrap(vars(name), scales = "free", strip.position = "bottom") +
  xlab(NULL) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.title = element_text(size = 10, face= "bold"))+
  theme(axis.text.x = element_text(color="black", size=8, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"))

##subset day 33
data_trt33 <- subset(data2, day == "33" )
str(data_trt33)
data_trt33$name <- as.factor(data_trt33$name)
levels(data_trt33$name) <- list("Bacteroides"="Bacteroides", "Elusimicrobium"="Elusimicrobium", "Neisseria"="Neisseria",
                                "Moraxella"="Moraxella", "Planctomycetaceae unclass"="Planctomycetaceae unclass") 

levels(data_trt33$treatment) <- list("A"="A", "GLN"="GLN", "NA"="NA_C")
data_trt33$day <- as.factor(data_trt33$day)
str(data_trt33)

## ggplot for all the species
ggplot(data = data_trt33, aes(x=name, y=ave.abundance,  fill=treatment)) +
  geom_bar(aes(x=name, y=ave.abundance), position= position_dodge(1), stat="identity",  alpha=0.7) +
  theme_test() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  geom_errorbar( aes(x=name, ymin=ave.abundance-sd, ymax=ave.abundance+sd), position = position_dodge(1), stat = "identity", size=0.4, width=0.5) +
  theme(axis.text.x=element_text(colour="black")) + 
  scale_y_continuous() +
  theme(axis.text.y=element_text(colour="black")) + 
  xlab("Genus") + ylab("Relative Abundance") + 
  labs(fill= "Treatment") +
  facet_wrap(vars(name), scales = "free", strip.position = "bottom") +
  xlab(NULL) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.title = element_text(size = 10, face= "bold"))+
  theme(axis.text.x = element_text(color="black", size=8, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"))
