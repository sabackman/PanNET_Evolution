library("ggplot2")
library("cowplot")
library("ggbreak")


# PanNET 1
samples1 <- rep(c("M2","P1c","M1b","P1a","M1a","P1b"),9)
genes1 <- c(rep("MEN1", 6),rep("ATRX", 6),rep("DAXX", 6),rep("TP53", 6),
            rep("PTEN", 6),rep("DEPDC5", 6),rep("TSC2", 6),rep("ARID1A", 6),rep("SETD2", 6))
mut1 <- c(rep("None", 6*9))
df1 <- data.frame(Sample=samples1,Gene=genes1, Mutation=mut1)

none_col="#EEEEEE"
  fs_col="#284f4d"
    indel_col="#796276"
      missense_col="#c371a2"
        stopgain_col="#cda6ae"
          splice_col="#d2d9b8"
            
          
          #PanNET2
          samples2 <- rep(c("P1b","P1a","M1"),9)
          genes2 <- c(rep("MEN1", 3),rep("ATRX", 3),rep("DAXX", 3),rep("TP53", 3),
                      rep("PTEN", 3),rep("DEPDC5", 3),rep("TSC2", 3),rep("ARID1A", 3),rep("SETD2", 3))
          mut2 <- c(rep("Indel", 3), c("Frameshift", "Frameshift", "Nonsense"), rep("None", 3),
                    rep("Missense", 3), rep("None", 5*3))
          df2 <- data.frame(Sample=samples2,Gene=genes2, Mutation=mut2)
          
          # PanNET 3
          samples3 <- rep(c("M4b", "P1a","P1b","M2","M4a","M3a","M1","M3b","P1c"),9)
          genes3 <- c(rep("MEN1", 9),rep("ATRX", 9),rep("DAXX", 9),rep("TP53", 9),
                      rep("PTEN", 9),rep("DEPDC5", 9),rep("TSC2", 9),rep("ARID1A", 9),rep("SETD2", 9))
          mut3 <- c(rep("Frameshift", 9), c("None", "Splice", "Splice", "None", "None", "None", "None", "None", "Splice"), c("Splice", "None", "None", "Splice", "Splice", "Splice", "Splice", "Splice", "None"), rep("None", 4*9), "Nonsense", "None", "None", rep("Nonsense", 4), "None", "None", "Frameshift", "None", "None", rep("Frameshift", 4), "None", "None")
          df3 <- data.frame(Sample=samples3,Gene=genes3, Mutation=mut3)
          
          ggplot(df3, aes(x=factor(Sample), y=factor(Gene, levels=c("MEN1", "ATRX", "DAXX", "PTEN", "DEPDC5", "SETD2", "TSC2", "TP53", "ARID1A")), fill=Mutation)) + geom_tile(color="#000000") + theme_classic() + scale_fill_manual(values=c(fs_col, none_col, splice_col, stopgain_col)) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position = "none")
          
          # PanNET 4
          samples4 <- rep(c("M1a","M1b","P1b","P1c","M1c","M3","P1a","M2b","M2a"),9)
          genes4 <- c(rep("MEN1", 9),rep("ATRX", 9),rep("DAXX", 9),rep("TP53", 9),
                      rep("PTEN", 9),rep("DEPDC5", 9),rep("TSC2", 9),rep("ARID1A", 9),rep("SETD2", 9))
          mut4 <- c(rep("Nonsense", 9), rep("None", 9), rep("Frameshift", 9), rep("None", 9), rep("Nonsense", 2), rep("None", 7), "None", "None", "Frameshift", "Frameshift", "None", "None", "None", "Frameshift", "None", rep("None", 3*9))
          df4 <- data.frame(Sample=samples4,Gene=genes4, Mutation=mut4)
          
          # PanNET 5
          samples5 <- rep(c("M1a", "M1b", "M2"), 9)
          genes5 <- c(rep("MEN1", 3),rep("ATRX", 3),rep("DAXX", 3),rep("TP53", 3),
                      rep("PTEN", 3),rep("DEPDC5", 3),rep("TSC2", 3),rep("ARID1A", 3),rep("SETD2", 3))
          mut5 <- c(rep("Indel", 3), rep("None", 3), rep("Missense", 3), rep("None", 3*3), rep("Missense", 3), rep("None", 3*2))
          
          # PanNET 6
          samples6 <- rep(c("M1", "M2"), 9)
          genes6 <- c(rep("MEN1", 2),rep("ATRX", 2),rep("DAXX", 2),rep("TP53", 2),
                      rep("PTEN", 2),rep("DEPDC5", 2),rep("TSC2", 2),rep("ARID1A", 2),rep("SETD2", 2))
          mut6 <- c(rep("None", 6*2), "Nonsense", "None", "Missense", "None", "None", "None")
          
          
          
          #
          # All samples
          #
          
          
          df_all <- data.frame(Sample=c(samples1, samples2, samples3, samples4, samples5, samples6), 
                               Gene=c(genes1, genes2, genes3, genes4, genes5, genes6), 
                               Mutation=c(mut1, mut2, mut3, mut4, mut5, mut6), 
                               Patient=c(rep("PanNET1", 54), rep("PanNET2", 27), rep("PanNET3", 81), rep("PanNET4", 81), rep("PanNET5", 27), rep("PanNET6", 18)))
          p2 <- ggplot(df_all, aes(x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=factor(Gene,levels=c("MEN1", "ATRX", "DAXX", "PTEN", "TSC2", "DEPDC5", "SETD2", "ARID1A", "TP53")), fill=Mutation)) + 
            geom_tile(color="#000000") + theme_classic()+
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(), 
                  strip.text.x = element_blank(), axis.ticks.x = element_blank(), 
                  axis.text.x = element_blank(), axis.line.x=element_blank()) + 
            scale_fill_manual(values=c(fs_col, indel_col,missense_col, none_col,splice_col , stopgain_col)) + 
            scale_x_discrete(expand = c(0, 0.5)) + 
            facet_grid(.~Patient, scale="free_x", space="free")
          
          
          Patients <- c(rep("PanNET4",9),rep("PanNET3",9),rep("PanNET2",3),rep("PanNET1",6), rep("PanNET5",3), rep("PanNET6",2))
          sample_name <- c("M1a","M1b","P1b","P1c","M1c","M3","P1a","M2b","M2a",
                           "M4b","P1a","P1b","M2","M4a","M3a","M1","M3b","P1c",
                           "P1b","P1a","M1",
                           "M2","P1c","M1b","P1a","M1a","P1b",
                           "M1a", "M1b", "M2",
                           "M1", "M2")
          mut_burden <- c(0.9610959762265735,0.9693046251338466,0.8369401615040659,1.0206086808043044,1.052075168282185,1.1505789551694638,0.9884581392508175,1.1714426044754498,1.21966841680568,
                          181.96590870605482,1.4416439643398602,1.7234742434895742,277.6992765871288,235.32794106297345,227.07995904635285,239.11794266886906,220.44121424259563,2.4126587246627222,
                          1.1878599022899963,1.158103550001131,1.1964105782350727,
                          2.7786276551119866,1.5784547794610806,2.1068865528667944,1.3014128788406092,2.485510483714772,1.3681081512122042,
                          
                          1.5376670075965875, 1.4472160071497295, 5.487360693776058, 
                          0.7236080035748648, 4.733602356718907)
          
          df_mut_burden <- data.frame(Sample=sample_name, Patient=Patients, MutBurden=mut_burden)
          
          p1 <- ggplot(df_mut_burden, 
                       aes(x=factor(Sample, levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=MutBurden)
          ) + geom_bar(stat="identity") + facet_grid(.~Patient, scale="free_x", space="free") + theme_classic() +
            scale_x_discrete(expand = c(0, 0.5)) + 
            theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
            labs(y="Mutation burden") + 
            scale_y_log10()
          
          timepoint <- c(rep("Baseline",5), "Followup", rep("Baseline", 3),
                         "Followup", rep("Baseline", 2), rep("Followup", 5), "Baseline",
                         "Baseline", "Baseline", "Followup",
                         "Followup", "Baseline", "Followup", "Baseline", "Followup", "Baseline", 
                         "Baseline", "Baseline", "Followup", 
                         "Baseline", "Followup")
          
          timepoint_df <- data.frame(Sample=sample_name, Patient=Patients, Timepoint=timepoint)
          
          p3 <- ggplot(timepoint_df, aes(x=factor(Sample, 
                                                  levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=1, fill=Timepoint)) +
            geom_tile(color="#000000") + facet_grid(.~Patient, scale="free_x", space="free") +
            theme_classic() + scale_x_discrete(expand = c(0, 0.5)) + 
            scale_fill_manual(values=c("#8da87a", "#dfdfc8")) +
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(), strip.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.direction="horizontal")
          
          ki67 <- c(1.8,1, 1.9, 1.7, 1, 3.3, NA, 1, 1, 
                    10, 1, 2.2, 60, 6.5, 75, 40, NA, 2.4,
                    4, 4.2, NA, 
                    4.3, 2, 3.1, 3.2, 4.5, 1.8, 16, 16, 20, 12, 30)
          ki67_df <- data.frame(Sample=sample_name, Patient=Patients, Ki67=ki67)
          
          p5 <- ggplot(ki67_df, aes(x=factor(Sample, 
                                             levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=1, fill=ki67)) +
            geom_tile(color="#000000") + facet_grid(.~Patient, scale="free_x", space="free") +
            theme_classic() + scale_x_discrete(expand = c(0, 0.5)) + 
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.line.x=element_blank(), strip.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.direction="horizontal") +
            scale_fill_gradient(low="#FFFFFF", high="#FF0000", breaks = c(20, 40, 60), labels=c("20%", "40%", "60%"))
          
          
          sig <- c("Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Present",
                   "Absent","Absent","Present","Present","Present","Present","Present","Absent","Absent","Absent",
                   "Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent","Absent",
                   "Absent","Absent")
          
          sig_df <- data.frame(Sample=sample_name, Patient=Patients, SBS11=sig)
          
          p_sig <- ggplot(sig_df, aes(x=factor(Sample, 
                                               levels=c("P1a", "P1b", "P1c", "M1", "M1a", "M1b", "M1c", "M2", "M2a", "M2b", "M3", "M3a", "M3b", "M4a", "M4b")), y=1, fill=SBS11)) +
            geom_tile(color="#000000") + facet_grid(.~Patient, scale="free_x", space="free") +
            theme_classic() + scale_x_discrete(expand = c(0, 0.5)) + 
            scale_fill_manual(values=c("#EEEEEE", "#FFA9A9")) +
            theme(axis.title.x = element_blank(),axis.title.y = element_blank(), axis.line.x=element_blank(), strip.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.direction="horizontal") 
          
          
          plot_grid(p1, p2, p5, p_sig, p3, nrow=5, align="v", axis="lr", rel_heights = c(1,2,0.3,0.3, 0.4))
          
          #
          #Statistical tests
          #
          #all_data <- data.frame(Sample=sample_name, Patients=Patients, TMB=mut_burden, Time=timepoint)
          
          #ggplot(all_data, aes(x=Time, y=log(TMB)))+ geom_point()
          