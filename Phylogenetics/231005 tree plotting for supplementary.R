library(ape)
library(ggtree)
library(ggplot2)
library(cowplot)

load("SNV_indel_trees.RData")
load("../PanNET_pairtree/Cloneplots.RData")

#Harmonize labels
colnames(data1)[colnames(data1) == "716c"] <- "P1c"
colnames(data1)[colnames(data1) == "5288"] <- "P1a"
colnames(data1)[colnames(data1) == "5292"] <- "P1b"
colnames(data1)[colnames(data1) == "10699"] <- "M1a"
colnames(data1)[colnames(data1) == "10701"] <- "M1b"
colnames(data1)[colnames(data1) == "10703"] <- "M2"
colnames(data2)[colnames(data2) == "2331.c"] <- "P1a"
colnames(data2)[colnames(data2) == "2331.d"] <- "P1b"
colnames(data2)[colnames(data2) == "A15447.6"] <- "M1"
colnames(data3)[colnames(data3) == "765.3"] <- "P1b"
colnames(data3)[colnames(data3) == "765.5"] <- "P1c"
colnames(data3)[colnames(data3) == "7370"] <- "M1"
colnames(data3)[colnames(data3) == "7372"] <- "M2"
colnames(data3)[colnames(data3) == "7379"] <- "M3a"
colnames(data3)[colnames(data3) == "7382"] <- "M3b"
colnames(data3)[colnames(data3) == "7385"] <- "M4a"
colnames(data3)[colnames(data3) == "7386"] <- "M4b"
colnames(data3)[colnames(data3) == "5621"] <- "P1a"
colnames(data4)[colnames(data4) == "2580"] <- "P1b"
colnames(data4)[colnames(data4) == "2584"] <- "M1a"
colnames(data4)[colnames(data4) == "2586"] <- "M1b"
colnames(data4)[colnames(data4) == "320.lever"] <- "M1c"
colnames(data4)[colnames(data4) == "320"] <- "P1a"
colnames(data4)[colnames(data4) == "6680"] <- "M3"
colnames(data4)[colnames(data4) == "2597"] <- "M2a"
colnames(data4)[colnames(data4) == "2600"] <- "M2b"
colnames(data4)[colnames(data4) == "320.10"] <- "P1c"

#data1$Germline <- NULL
#data2$Germline <- NULL
#data3$Germline <- NULL
#data4$Germline <- NULL

dist1 <- dist(t(data1))
dist2 <- dist(t(data2))
dist3 <- dist(t(data3))
dist4 <- dist(t(data4))
evo1 <- fastme.bal(dist1)
evo2 <- fastme.bal(dist2)
evo3 <- fastme.bal(dist3)
evo4 <- fastme.bal(dist4)

evo1 <- root(evo1, outgroup="Germline", resolve.root=TRUE)
evo2 <- root(evo2, outgroup="Germline", resolve.root=TRUE)
evo3 <- root(evo3, outgroup="Germline", resolve.root=TRUE)
evo4 <- root(evo4, outgroup="Germline", resolve.root=TRUE)

#ggtree(evo1) + geom_tiplab()
#ggtree(evo2) + geom_tiplab()
#ggtree(evo3) + geom_tiplab()
#ggtree(evo4) + geom_tiplab()


tree1 <- ggtree(evo1) + geom_tiplab() + 
  labs(title = "PanNET 1") + 
  theme(plot.title = element_text(hjust = 0.5))

ann_data4 <- data.frame(c(18, 13), c("PTEN subclonal", "DEPDC5"))
timedata4 <- data.frame(c(16,14,7), c("Baseline", "Baseline", "Follow-up"))
colnames(timedata4) <- c("node", "Timepoint")
colnames(ann_data4) <- c("node", "Variant")
tree4 <- ggtree(evo4) + geom_highlight(data=ann_data4, mapping=aes(node=node, fill=Variant), type="roundrect") + 
  geom_tiplab() + scale_fill_manual(values=c("#FFA69E", "#AED9E0")) +
  labs(title = "PanNET 4") + 
  theme(plot.title = element_text(hjust = 0.5))

ann_data2 <- data.frame(c(3, 6), c("Ser1913*", "Lys372fs"))
colnames(ann_data2) <- c("node", "ATRX")
tree2 <- ggtree(evo2) + geom_highlight(data=ann_data2, mapping=aes(node=node, fill=ATRX), type="roundrect") + geom_tiplab() + 
  scale_fill_manual(values=c("#F4E3B2", "#D3D5D7")) +
  labs(title = "PanNET 2") + 
  theme(plot.title = element_text(hjust = 0.5))

ann_data3 <- data.frame(c(1,2,9,14), c( "ATRX\nc.4121-2A>G",  "ATRX\nc.4121-2A>G", "ATRX\nc.4121-2A>G", "DAXX\nc.208-1G>A*"))
colnames(ann_data3) <- c("node", "Mutation")
tree3 <- ggtree(evo3) + geom_highlight(data=ann_data3, mapping=aes(node=node, fill=Mutation), type="roundrect") + geom_tiplab() + 
  scale_fill_manual(values=c("#F4E3B2", "#D3D5D7")) +
  labs(title = "PanNET 3") + 
  theme(plot.title = element_text(hjust = 0.5))


empty_plot <- ggplot() + theme_void()

ggsave("Tree1.tiff", tree1)
ggsave("Tree2.tiff", tree2)
ggsave("Tree3.tiff", tree3)
ggsave("Tree4.tiff", tree4)

