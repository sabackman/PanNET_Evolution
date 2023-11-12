library(ggplot2)
library(cowplot)

data <- read.table("data.txt", sep="\t", header=T)
data$TMB_log2 <- log2(data$TMB)
data[data$HighGradeProgression==1,]$HighGradeProgression <- "High grade progression"
data[data$HighGradeProgression==0,]$HighGradeProgression <- "No high grade progression"

#
# Split data for stacked barplots
#

n <- length(data$Patient_ID)*3
sample <- c(rep("", n))
type <- c(rep("", n))
TMB <- c(rep(0, n))
HGP <- c(rep("", n))

i <- 1
for (p_id in data$Patient_ID) {
  sample[i] <- p_id
  type[i] <- "Alkylating"
  HGP[i] <- data[data$Patient_ID==p_id,]$HighGradeProgression
  if (!is.na(data[data$Patient_ID==p_id,]$Proportion_alkylating)) {
    TMB[i] <- data[data$Patient_ID==p_id,]$Proportion_alkylating*data[data$Patient_ID==p_id,]$TMB
  }
  
  i <- i+1
  sample[i] <- p_id
  type[i] <- "Not alkylating"
  HGP[i] <- data[data$Patient_ID==p_id,]$HighGradeProgression
  if (!is.na(data[data$Patient_ID==p_id,]$Proportion_alkylating)) {
    TMB[i] <- (1-data[data$Patient_ID==p_id,]$Proportion_alkylating)*data[data$Patient_ID==p_id,]$TMB
  }
  
  i <- i+1
  sample[i] <- p_id
  type[i] <- "NA"
  HGP[i] <- data[data$Patient_ID==p_id,]$HighGradeProgression
  if (is.na(data[data$Patient_ID==p_id,]$Proportion_alkylating)) {
    TMB[i] <- data[data$Patient_ID==p_id,]$TMB
  }
  i <- i+1
}

processed_data <- data.frame(sample, HGP, type, TMB)


#
# Prepare MMR gene status for plot
#
data$MSH6 <- c(rep(0, n/3))
data$MSH2 <- c(rep(0, n/3))
data$MLH3 <- c(rep(0, n/3))
for (p in data$Patient_ID) {
  if (grepl("MSH6", data[data$Patient_ID==p,]$MMR_Gene_variants, fixed=T)) {
    data[data$Patient_ID==p,]$MSH6 <- 1
  }
  if (grepl("MSH2", data[data$Patient_ID==p,]$MMR_Gene_variants, fixed=T)) {
    data[data$Patient_ID==p,]$MSH2 <- 1
  }
  if (grepl("MLH3", data[data$Patient_ID==p,]$MMR_Gene_variants, fixed=T)) {
    data[data$Patient_ID==p,]$MLH3 <- 1
  }
}

data$MSH6 <- factor(data$MSH6)
data$MSH2 <- factor(data$MSH2)
data$MLH3 <- factor(data$MLH3)

mut <- c(data$MSH6, data$MSH2, data$MLH3)
gene <- c(rep("MSH6", n/3), rep("MSH2", n/3), rep("MLH3", n/3))
p_id <- rep(data$Patient_ID, 3)
HGP <- rep(data$HighGradeProgression ,3)
mut_df <- data.frame(p_id, HGP, gene, mut)

p1 <- ggplot(processed_data, aes(fill=type, y=TMB, x=factor(sample, levels=rev(data[order(data$TMB),]$Patient_ID)))) + 
  geom_bar(position="stack", stat="identity") + facet_grid(~HGP, space="free", scales="free", drop=TRUE) + 
  theme_classic() + 
  scale_fill_manual(values=c("#4b42f5", "#a3a3a3", "#f5bf42"), name="Mutation origin", labels=c("Alkylation related", "Unknown", "Other")) +
  theme(axis.title.x = element_blank(), axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  geom_hline(yintercept=50, color="#000000") 


p_mut <- ggplot(mut_df, aes(x=factor(p_id, levels=rev(data[order(data$TMB),]$Patient_ID)), fill=mut, y=gene)) + geom_tile() + 
  facet_grid(~HGP, space="free", scales="free") + 
  scale_fill_manual(values=c("#FFFFFF", "#000000"), name="Mutation", labels=c("No", "Yes")) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.text.y=element_text(face="italic"), strip.text.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), legend.direction="vertical") +
  xlab("Sample")

plot_grid(p1, p_mut, nrow=2, align="v", axis="lr", rel_heights=c(1.5, 0.5)) 
