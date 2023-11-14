library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(clusterProfiler)
library("org.Hs.eg.db", character.only=TRUE)

# All_samples
samples <- read.table("samples.txt", sep="\t", header=T)
samples$type <- relevel(factor(samples$type, levels=c("P", "M")), ref="P")
samples$patient <- factor(samples$patient)

files <- file.path("FAT_star_salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("FAT_star_salmon", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ patient + type)
dds <- DESeq(ddsTxi)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]

EnhancedVolcano(res, rownames(res), x="log2FoldChange", y="padj")

res_fc <- res$stat
names(res_fc) <- rownames(res)
res_fc <- na.omit(res_fc)
res_fc <- sort(res_fc, decreasing=TRUE)

res_fc <- data.frame(res_fc) %>% mutate(rank = rank(res_fc,  ties.method = "random")) %>% arrange(desc(rank))

res_fc$rankfixed <- res_fc$rank - (length(res_fc$res_fc)-sum(ifelse(res_fc$res_fc > 0, 1, 0))) -1

gene_list <- res_fc$rankfixed
names(gene_list) <- rownames(res_fc)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL",
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             eps = 0,
             pAdjustMethod = "BH",
             nPermSimple = 100000)
resnona <- na.omit(res)
up_enrich <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange> 1.5,]), org.Hs.eg.db, keyType="SYMBOL", ont="BP")
down_enrich <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -1.5,]), org.Hs.eg.db, keyType="SYMBOL", ont="BP")










#PanNET 3
samples <- read.table("PanNET3_samples.txt", sep="\t", header=T)

files <- file.path("FAT_star_salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("FAT_star_salmon", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

samples$type <- relevel(factor(samples$type, levels=c("P", "M")), ref="P")

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ type)
dds <- DESeq(ddsTxi)
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="PanNET3_results.csv")
EnhancedVolcano(res, rownames(res), x="log2FoldChange", y="padj")

#GSEA
res_fc <- res$stat
names(res_fc) <- rownames(res)
res_fc <- na.omit(res_fc)
res_fc <- sort(res_fc, decreasing=TRUE)

res_fc <- data.frame(res_fc) %>% mutate(rank = rank(res_fc,  ties.method = "random")) %>% arrange(desc(rank))

res_fc$rankfixed <- res_fc$rank - (length(res_fc$res_fc)-sum(ifelse(res_fc$res_fc > 0, 1, 0))) -1

gene_list <- res_fc$rankfixed
names(gene_list) <- rownames(res_fc)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL",
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             eps = 0,
             pAdjustMethod = "BH",
             nPermSimple = 100000)

resnona <- na.omit(res)
up_enrich <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange> 1.5,]), org.Hs.eg.db, keyType="SYMBOL", ont="BP")
down_enrich <- enrichGO(rownames(resnona[resnona$padj<0.05 & resnona$log2FoldChange< -1.5,]), org.Hs.eg.db, keyType="SYMBOL", ont="BP")

up_enrich <- mutate(up_enrich, q=-log(p.adjust, base=10))
write.table(up_enrich, "PanNET3_GO_up.txt")
down_enrich <- mutate(down_enrich, q=-log(p.adjust, base=10))
write.table(down_enrich, "PanNET3_GO_down.txt")

barplot(up_enrich, x="q", showCategory = 12) + 
  theme_classic() + xlab("-log10 Adjusted p-value")
barplot(down_enrich, x="q", showCategory = 12) + theme_classic() +
  xlab("-log10 Adjusted p-value")
