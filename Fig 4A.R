library("ggplot2")
library("cowplot")
library("ggbreak")
library("MutationalPatterns")

# PanNET 1
samples1 <- rep(c("M2","P1c","M1b","P1a","M1a","P1b"),9)
df1 <- data.frame(Sample=samples1)

#PanNET2
samples2 <- rep(c("P1b","P1a","M1"),9)
df2 <- data.frame(Sample=samples2)

# PanNET 3
samples3 <- rep(c("M4b", "P1a","P1b","M2","M4a","M3a","M1","M3b","P1c"),9)
df3 <- data.frame(Sample=samples3)

# PanNET 4
samples4 <- rep(c("M1a","M1b","P1b","P1c","M1c","M3","P1a","M2b","M2a"),9)
df4 <- data.frame(Sample=samples4)

#
# All samples
#


Patients <- c(rep("PanNET4",9),rep("PanNET3",9),rep("PanNET2",3),rep("PanNET1",6))
sample_name <- c("M1a","M1b","P1b","P1c","M1c","M3","P1a","M2b","M2a",
                 "M4b","P1a","P1b","M2","M4a","M3a","M1","M3b","P1c",
                 "P1b","P1a","M1",
                 "M2","P1c","M1b","P1a","M1a","P1b")
mut_burden <- c(0.95,0.96,0.83,1.02,1.04,1.14,0.98,1.16,1.21,
                181.61,1.43,1.71,277.08,234.88,226.60,238.59,219.98,2.40,
                1.18,1.15,1.19,
                2.77,1.57,2.10,1.29,2.47,1.36)

df_mut_burden <- data.frame(Sample=sample_name, Patient=Patients, MutBurden=mut_burden)

p1 <- ggplot(df_mut_burden, 
             aes(x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=MutBurden)
) + geom_bar(stat="identity", fill="#326fa8") + facet_grid(.~Patient, scale="free_x", space="free") + theme_classic() +
  scale_x_discrete(expand = c(0, 0.5)) + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  labs(y="Mutation burden") + 
  scale_y_log10()

#
#MSH6 plot
#
muts <- c("p.Ala1226Val","p.Asp1255Asn","p.Gln1155*","p.Gly409Glu","p.Pro1097Leu","p.Pro386Ser","p.Pro399Leu","p.Thr1008Ile")

empty_VAF <- c(0,0,0,0,0,0,0,0)
S7370 <- c(0.272,0.159,0.181,0,0,0,0.258,0)
S7372 <- c(0.245,0,0.333,0,0,0.238,0.256,0)
S7379 <- c(0.471,0,0.409,0,0,0.571,0.568,0)
S7382 <- c(0.410,0,0.395,0,0,0.393,0.459,0)
S7385 <- c(0,0,0,0,0,0,0,0.065)
S7386 <- c(0,0,0,0.230,0.300,0,0,0)

d = data.frame("Mutation"=rep(muts, 6),
               "Sample"=c(rep("M1",8),rep("M2",8),rep("M3a",8),rep("M3b",8),rep("M4a",8),rep("M4b",8)),
               "VAF"= c(S7370, S7372, S7379, S7382, S7385, S7386))

samples1_VAF <- rep(c("M2","P1c","M1b","P1a","M1a","P1b"),8)
samples2_VAF <- rep(c("P1b","P1a","M1"),8)
samples4_VAF <- rep(c("M1a","M1b","P1b","P1c","M1c","M3","P1a","M2b","M2a"),8)


d = data.frame("Mutation"=rep(muts, 27),
               "Patient"=c(rep("PanNET1",6*8), rep("PanNET2",3*8), rep("PanNET3",9*8), rep("PanNET4",9*8)), 
               "Sample"=c(samples1_VAF, samples2_VAF, rep("P1a", 8), rep("P1b", 8), rep("P1c", 8), rep("M1",8),rep("M2",8),rep("M3a",8),rep("M3b",8),rep("M4a",8),rep("M4b",8), samples4_VAF),
               "VAF"= c(rep(empty_VAF, 6), rep(empty_VAF, 3), rep(empty_VAF, 3), S7370, S7372, S7379, S7382, S7385, S7386, rep(empty_VAF,9)))

p_mut <- ggplot(d, aes(x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=factor(Mutation, levels=c("p.Ala1226Val","p.Pro399Leu","p.Gln1155*","p.Pro386Ser","p.Asp1255Asn","p.Thr1008Ile","p.Gly409Glu","p.Pro1097Leu")), fill=VAF)) + 
  geom_tile() + scale_fill_gradient(low="#FFFFFF", high="#FF0000") +
  facet_grid(.~Patient, scale="free_x", space="free") +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x=element_blank())


#
# Signature plot
#
load("2_nmf_res.RData")
Sig_A <- data.frame(t(nmf_res$contribution))$Signature.A[2:28]
Sig_B <- data.frame(t(nmf_res$contribution))$Signature.B[2:28]
samples <- rownames(data.frame(t(nmf_res$contribution)))[2:28]
Patients <- c(rep("PanNET1",6), rep("PanNET2",3), rep("PanNET3",9), rep("PanNET4",9))
samples <- c("M1a","M1b","M2","P1a","P1b","P1c",
             "P1a","P1b","M1",
             "P1a","M1","M2","M3a","M3b","M4a","M4b","P1b","P1c",
             "P1b","M1a","M1b","M2a","M2b","P1a","P1c","M1c","M3")

d = data.frame("Sample"=samples,
               "Patient"=Patients,
               "Signature A"=Sig_A,
               "Signature B"=Sig_B)

d$Contribution <- d$Signature.A/(d$Signature.B+d$Signature.A)

p_sig <- ggplot(d, 
       aes(y=Contribution, x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")))) + 
  geom_bar(stat="identity", fill="#44a832") +
  facet_grid(.~Patient, scale="free_x", space="free") +
  theme_classic() +
  theme(axis.title.x = element_blank(),strip.text.x=element_blank()) + 
  labs(y="Relative contribution\nof Signature A")

#
# MSI plot
#

pat <- c("PanNET1","PanNET1","PanNET1","PanNET1","PanNET1","PanNET1","PanNET2","PanNET2","PanNET2","PanNET3","PanNET3","PanNET3","PanNET3","PanNET3","PanNET3","PanNET3","PanNET3","PanNET3","PanNET4","PanNET4","PanNET4","PanNET4","PanNET4","PanNET4","PanNET4","PanNET4","PanNET4")
sample <- c("M1b","P1a","M1a","M2","P1b","P1c","P1a","P1b","M1","M3a","P1a","M3b","P1c","M4b","P1b","M2","M1","M4a","M2b","M2a","M1c","P1a","M1a","P1c","M1b","P1b","M3")
MSI <- c(0.05,0.05,0.04,0.02,0,0,0.26,0.2,0.09,2.81,2.07,1.83,1.4,1.09,0.93,0.52,0.19,0,2.87,2.67,2.02,1.69,1.56,1.54,1.28,0.62,0.54)

msidata <- data.frame(Patient=pat, MSIsensor = MSI, Sample <- sample)

plot_msi <- 
  ggplot(msidata, aes(x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=MSIsensor)) + 
    geom_bar(stat="identity", fill="#a87d32") + 
    facet_grid(.~Patient, scale="free_x", space="free") +
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x=element_blank()) + 
  scale_y_continuous(limits=c(0,4)) + 
  geom_hline(yintercept=3.5, color="#777777") 

plot_grid(p1, p_sig, p_mut, plot_msi, nrow=4, align="v", axis="lr", rel_heights = c(1,1,1))

