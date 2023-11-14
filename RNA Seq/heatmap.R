library(DESeq2)
library(tximport)
library(gplots)

variation_measure <- function (x) {
  if (sum(x) == 0) {
    return(0)
  }
  return(var(x)/mean(x))
}

# All_samples
samples <- read.table("samples.txt", sep="\t", header=T)

files <- file.path("FAT_star_salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read.table(file.path("FAT_star_salmon", "salmon_tx2gene.tsv"))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ patient + type)
dds <- DESeq(ddsTxi)

counts_normalized <- as.data.frame(counts(dds, normalized=TRUE))
variation <- apply(counts_normalized, 1, variation_measure)
counts_top <- counts_normalized[variation>quantile(variation, c(0.95)),]

counts_top_scaled <- t(scale(t(counts_top)))

#Fix sample names for publication
#> colnames(counts_top_scaled)
#[1] "SJ-2336-10699"     "SJ-2336-10701"     "SJ-2336-10703"     "SJ-2336-2331-c"   
#[5] "SJ-2336-2331-d"    "SJ-2336-2580"      "SJ-2336-2584"      "SJ-2336-2586"     
#[9] "SJ-2336-2597"      "SJ-2336-320"       "SJ-2336-320-10"    "SJ-2336-320-lever"
#[13] "SJ-2336-5288"      "SJ-2336-5292"      "SJ-2336-5621"      "SJ-2336-6680"     
#[17] "SJ-2336-716c"      "SJ-2336-7370"      "SJ-2336-7372"      "SJ-2336-7379"     
#[21] "SJ-2336-7382"      "SJ-2336-7385"      "SJ-2336-7386"      "SJ-2336-765-3"    
#[25] "SJ-2336-765-5"     "SJ-2336-A15447-06"
colnames(counts_top_scaled) <- c("M1a", "M1b", "M2", "P1a",
                                 "P1b", "P1b", "M1a", "M1b",
                                 "M2a", "P1a", "P1c", "M1c",
                                 "P1a", "P1b", "P1a", "M3",
                                 "P1c", "M1", "M2", "M3a",
                                 "M3b", "M4a", "M4b", "P1b",
                                 "P1c", "M1")
rsc <- samples$patient
rsc[rsc==1] <- "#FF0000"
rsc[rsc==2] <- "#0000FF"
rsc[rsc==3] <- "#00FF00"
rsc[rsc==4] <- "#FFFF00"
heatmap.2(t(counts_top_scaled), trace="none", 
          density.info ="none", col=colorRampPalette(c('green',"black",'red'))(32),
          labCol=NA, RowSideColors = rsc)

