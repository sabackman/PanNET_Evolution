library(minfi)
library(conumee)
library(gplots)

baseDir <- "SD-2151_191218_IDAT"

targets <- read.metharray.sheet(baseDir, pattern="pnet$")
RGset <- read.metharray.exp(targets = targets)

methylSet <- preprocessIllumina(RGset)
ratioSet <- ratioConvert(methylSet, what="both")
gset <- mapToGenome(ratioSet)

beta <- getBeta(gset)
write.table(beta, "beta_values_pnet.txt", sep="\t")

#
# Extract MGMT
# promoter methylation
#
annotation <- getAnnotation(gset)

annotation <- annotation[annotation$UCSC_RefGene_Name %in% c("MGMT", "MGMT;MGMT"),]

annotation_small <- data.frame(annotation$chr)
annotation_small$pos <- annotation$pos
annotation_small$Name <- annotation$Name
annotation_small$Relation_to_Island <- annotation$Relation_to_Island

#Full MGMT region
MGMT_meth <- beta[rownames(beta) %in% rownames(annotation),]
colnames(MGMT_meth) <- targets$Sample_Name
heatmap(MGMT_meth)

#CpG island
MGMT_meth <- beta[rownames(beta) %in% rownames(annotation[annotation$Relation_to_Island=="Island",]),]
colnames(MGMT_meth) <- targets$Sample_Name
heatmap(MGMT_meth, Rowv=NA, Colv=NA)

#
#
# Phylogenetics 
#
#
library(minfi)
library(ape)
library(gplots)

baseDir <- "SD-2151_191218_IDAT"

targets <- read.metharray.sheet(baseDir, pattern="pnet$")
targets <- targets[targets$patient!=7,] #Remove single PNET from MSH project

for (patient in unique(targets$patient)) {

  targets_current <- targets[targets$patient==patient,]
  
  RGset <- read.metharray.exp(targets = targets_current)
  
  methylSet <- preprocessIllumina(RGset)
  ratioSet <- ratioConvert(methylSet, what="both")
  gset <- mapToGenome(ratioSet)
  
  beta <- getBeta(gset)
  for (name in colnames(beta)) {
    colnames(beta)[colnames(beta)==name] <- paste(targets_current[targets_current$Array==strsplit(name, "_")[[1]][2] & targets_current$Slide==strsplit(name, "_")[[1]][1],]$Sample_Name, 
                                                  substring(targets_current[targets_current$Array==strsplit(name, "_")[[1]][2] & targets_current$Slide==strsplit(name, "_")[[1]][1],]$Type[1], 1, 1), 
                                                  sep="_")
  }
  
  #Remove probes with missing values
  beta <- beta[complete.cases(beta),]
  
  #Find 1% of probes with largest interval
  diff <- apply(beta, 1, function(x) max(x)-min(x))
  dm <- dist(t(beta[diff>quantile(diff, c(0.99)),]))
  save(dm, file=paste(patient, "_methdist.RData", sep=""))
  evo <- fastme.bal(dm)
  save(evo, file=paste(patient, "_tree.RData", sep=""))
}

#
#
# Supplementary figures
#
#

met_color <- "#FF6666"
primary_color <- "#6666FF"
tiplabs = TRUE

#
# Set up colors for PanNET 3
#
load("PanNET_03_tree.RData")
colors <- c(rep("black", 15))
colors[13] <- met_color
colors[12] <- met_color
colors[11] <- met_color
colors[10] <- met_color
colors[8] <- met_color
colors[6] <- met_color

colors[1] <- primary_color
colors[15] <- primary_color
colors[14] <- primary_color

tips <- c("black", "grey", "grey", "grey", "grey", "grey", "grey", "black", "black")
evo$tip.label <- c("P1a","M1","M2","M3a","M3b","M4b","M4a","P1b","P1c")

plot(evo, "unrooted", edge.color=colors, tip.color=tips, edge.width=3, main="PanNET 3",
     show.tip.label=tiplabs, rotate.tree=110)


#
# Set up colors for PanNET 1
#
load("PanNET_01_tree.RData")
colors <- c(rep("black", 9))
colors[6] <- met_color
colors[8] <- met_color
colors[9] <- met_color

colors[4] <- primary_color
colors[3] <- primary_color
colors[1] <- primary_color

tips <- c("black", "black", "black", "grey", "grey", "grey")

evo$tip.label <- c("P1a","P1b","P1c","M1a","M1b","M2")

plot(evo, "unrooted", edge.color=colors, tip.color=tips, edge.width=3, main="PanNET 1",
     show.tip.label=tiplabs)


#
# Set up colors for PanNET 2
#
load("PanNET_02_tree.RData")
colors <- c(rep(met_color, 3))

tips <- c("black", "black", "grey")
evo$tip.label <- c("P1a","P1b","M1")


plot(evo, "unrooted", edge.color=colors, tip.color=tips, edge.width=3, main="PanNET 2",
     show.tip.label=tiplabs, rotate.tree=-28)

#
# Set up colors for PanNET 4
#
load("PanNET_04_tree.RData")
colors <- c(rep("black", 15))

colors[15] <- met_color
colors[14] <- met_color
colors[13] <- met_color
colors[11] <- met_color
colors[4] <- met_color
colors[3] <- met_color
colors[12] <- primary_color
colors[9] <- primary_color
colors[1] <- primary_color

tips <- c("black", "black", "black", "black", "black", "black", "black", "grey", "black")
evo$tip.label <- c("P1a", "M1a", "M1b", "P1b", "M2b","P1c","M2a",
                   "M3", "M1c")

plot(evo, "unrooted", edge.color=colors, tip.color=tips, edge.width=3, main="PanNET 4",
     show.tip.label=tiplabs, rotate.tree=-36)



#
# Differential methylation for PanNET3
#
library(minfi)
library(gplots)
library(ggplot2)

baseDir <- "SD-2151_191218_IDAT"

targets <- read.metharray.sheet(baseDir, pattern="pnet$")
targets <- targets[targets$patient=="PanNET_03",] 

RGset <- read.metharray.exp(targets = targets)

methylSet <- preprocessIllumina(RGset)
ratioSet <- ratioConvert(methylSet, what="both")
gset <- mapToGenome(ratioSet)

beta <- getBeta(gset)
for (name in colnames(beta)) {
  colnames(beta)[colnames(beta)==name] <- paste(targets[targets$Array==strsplit(name, "_")[[1]][2] & targets$Slide==strsplit(name, "_")[[1]][1],]$Sample_Name, 
                                                substring(targets[targets$Array==strsplit(name, "_")[[1]][2] & targets$Slide==strsplit(name, "_")[[1]][1],]$Type[1], 1, 1), 
                                                sep="_")
}

#Remove probes with missing values
beta <- beta[complete.cases(beta),]
#
# Check difference between primary and met
#
betadiff <- apply(beta[,targets$Type!="Primary"], 1, mean)-apply(beta[,targets$Type=="Primary"], 1, mean)
ann <- getAnnotation(gset, lociNames=names(betadiff))


d <- data.frame(delta=abs(betadiff), type=ifelse(ann$Relation_to_Island=="Island", "Island", "Other"))
ttest <- t.test(d[d$type=="Island",]$delta, d[d$type!="Island",]$delta)


#
# Supplementary figure
#
ggplot(d, aes(x=delta, color=type)) + geom_density() + theme_classic() +
  xlab("Abs(delta)") + ylab("Density") + ggtitle("Mean Delta between PanNET3 \nprimary and metastatic samples") +
  annotate("text", x=0.6, y=20, label="p-value < 2.2e-16") + 
  labs(color='CpG type') 

