# Set working directory and load packages 
setwd("~/Desktop/MIARNAseq code and plots")
library(lattice)
library(DESeq2) 
library(pheatmap) 
library(GSA) 
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(gprofiler2)
library(pals)
library(EnhancedVolcano)
library(stringr)

#we create a function to define our genes as not changed, or up or downregulated, this is with padj
categorize.deseq.df <- function(df, thresh = 0.1, log2fold = 0.0, treat = 'Auxin') {
  df.activated = data.frame(matrix(nrow = 0, ncol = 0))
  df.repressed = data.frame(matrix(nrow = 0, ncol = 0))
  if (nrow(df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold,]) != 0) {
    df.activated = df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold,]
    df.activated$response = paste(treat, 'Activated')
  }
  if (nrow(df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold,]) != 0) {
    df.repressed = df[df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold,]
    df.repressed$response = paste(treat, 'Repressed')
  }
  df.unchanged = df[df$padj > 0.5 & !is.na(df$padj) & abs(df$log2FoldChange) < 0.25,]
  df.unchanged$response = paste(treat, 'Unchanged')
  df.dregs = df[!(df$padj < thresh & !is.na(df$padj) & df$log2FoldChange > log2fold) &
                  !(df$padj < thresh & !is.na(df$padj) & df$log2FoldChange < -log2fold) &
                  !(df$padj > 0.5 & !is.na(df$padj) &
                      abs(df$log2FoldChange) < 0.25), ]
  df.dregs$response = paste(treat, 'All Other Genes')
  df.effects.lattice =
    rbind(df.activated,
          df.unchanged,
          df.repressed,
          df.dregs)
  df.effects.lattice$response = factor(df.effects.lattice$response)
  df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'Unchanged'))
  df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'All Other Genes'))
  return(df.effects.lattice)
}

#this is with pval 
categorize.deseq.df2 <- function(df, thresh = 0.05, log2fold = 0.0, treat = 'Auxin') {
  df.activated = data.frame(matrix(nrow = 0, ncol = 0))
  df.repressed = data.frame(matrix(nrow = 0, ncol = 0))
  if (nrow(df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold,]) != 0) {
    df.activated = df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold,]
    df.activated$response = paste(treat, 'Activated')
  }
  if (nrow(df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold,]) != 0) {
    df.repressed = df[df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold,]
    df.repressed$response = paste(treat, 'Repressed')
  }
  df.unchanged = df[df$pval > 0.5 & !is.na(df$pval) & abs(df$log2FoldChange) < 0.25,]
  df.unchanged$response = paste(treat, 'Unchanged')
  df.dregs = df[!(df$pval < thresh & !is.na(df$pval) & df$log2FoldChange > log2fold) &
                  !(df$pval < thresh & !is.na(df$pval) & df$log2FoldChange < -log2fold) &
                  !(df$pval > 0.5 & !is.na(df$pval) &
                      abs(df$log2FoldChange) < 0.25), ]
  df.dregs$response = paste(treat, 'All Other Genes')
  df.effects.lattice2 =
    rbind(df.activated,
          df.unchanged,
          df.repressed,
          df.dregs)
  df.effects.lattice2$response = factor(df.effects.lattice2$response)
  df.effects.lattice2$response = relevel(df.effects.lattice2$response, ref = paste(treat, 'Unchanged'))
  df.effects.lattice2$response = relevel(df.effects.lattice2$response, ref = paste(treat, 'All Other Genes'))
  return(df.effects.lattice2)
}

#PCA plot function
plotPCAlattice <- function(df, file = 'PCA_lattice.pdf') {
  perVar = round(100 * attr(df, "percentVar"))
  df = data.frame(cbind(df, sapply(strsplit(as.character(df$name), '_rep'), '[', 1)))
  colnames(df) = c(colnames(df)[1:(ncol(df)-1)], 'unique_condition')
  print(df)
  color.x = substring(rainbow(length(unique(df$unique_condition))), 1,7)
  df$color = NA
  df$alpha.x = NA
  df$alpha.y = NA
  df$colpal = NA
  for (i in 1:length(unique(df$unique_condition))) {
    df[df$unique_condition == unique(df$unique_condition)[[i]],]$color = color.x[i]
    reps_col<- df[df$unique_condition == unique(df$unique_condition)[[i]],]
    replicates.x = nrow(reps_col)
    alx <- rev(seq(0.2, 1, length.out = replicates.x))
    for(rep in 1:replicates.x) {
      na <- reps_col[rep, ]$name
      df[df$name == na, ]$alpha.x = alx[rep]
      aly = as.hexmode(round(alx * 255))
      df[df$name == na, ]$alpha.y = aly[rep]
      cp = paste0(color.x[i], aly)
      df[df$name == na, ]$colpal = cp[rep]
    }
  }
  colpal = df$colpal
  df$name = gsub('_', ' ', df$name)
  df$name <- factor(df$name, levels=df$name, order=TRUE)
  pdf(file, width=6, height=6, useDingbats=FALSE)
  print(xyplot(PC2 ~ PC1, groups = name, data=df,
               xlab = paste('PC1: ', perVar[1], '% variance', sep = ''),
               ylab = paste('PC2: ', perVar[2], '% variance', sep = ''),
               par.settings = list(superpose.symbol = list(pch = c(20), col=colpal)),
               pch = 20, cex = 1.7,
               aspect = 1,
               auto.key = TRUE,
               col = colpal))
  dev.off()
}


#now we are loading the data in and creating the dataframe
x <- read.table("MG_PolyIC_E12_F1.gene.counts.txt",header=FALSE)

all.counts <- data.frame(row.names = x$V1)

samples <- c("MG_PolyIC_E12_F1", "MG_PolyIC_E12_F2","MG_PolyIC_E12_F3","MG_PolyIC_E12_F4",
             "MG_PolyIC_E12_M1","MG_PolyIC_E12_M2","MG_PolyIC_E12_M3","MG_PolyIC_E12_M4",
             "MG_PolyIC_E14_F1","MG_PolyIC_E14_F2","MG_PolyIC_E14_F3","MG_PolyIC_E14_F4",
             'MG_PolyIC_E14_M1','MG_PolyIC_E14_M2','MG_PolyIC_E14_M3','MG_PolyIC_E14_M4',
             "MG_PolyIC_P5_F1", "MG_PolyIC_P5_F2","MG_PolyIC_P5_F3","MG_PolyIC_P5_F4",
             "MG_PolyIC_P5_M1", "MG_PolyIC_P5_M2","MG_PolyIC_P5_M3","MG_PolyIC_P5_M4",
             "MG_PolyIC_P90_F1", "MG_PolyIC_P90_F2","MG_PolyIC_P90_F3","MG_PolyIC_P90_F4",
             "MG_PolyIC_P90_M1", "MG_PolyIC_P90_M2","MG_PolyIC_P90_M3","MG_PolyIC_P90_M4",
             "MG_Saline_E12_F1", "MG_Saline_E12_F2","MG_Saline_E12_F3","MG_Saline_E12_F4",
             "MG_Saline_E12_M1","MG_Saline_E12_M2","MG_Saline_E12_M3","MG_Saline_E12_M4",
             "MG_Saline_E14_F1","MG_Saline_E14_F2","MG_Saline_E14_F3","MG_Saline_E14_F4",
             'MG_Saline_E14_M1','MG_Saline_E14_M2','MG_Saline_E14_M3','MG_Saline_E14_M4',
             "MG_Saline_P5_F1", "MG_Saline_P5_F2","MG_Saline_P5_F3","MG_Saline_P5_F4",
             "MG_Saline_P5_M1", "MG_Saline_P5_M2","MG_Saline_P5_M3","MG_Saline_P5_M4",
             "MG_Saline_P90_F1", "MG_Saline_P90_F2","MG_Saline_P90_F3","MG_Saline_P90_F4",
             "MG_Saline_P90_M1", "MG_Saline_P90_M2","MG_Saline_P90_M3","MG_Saline_P90_M4")




for (i in samples){
  j<-read.table(print(paste0(i,".gene.counts.txt")))
  all.counts  = cbind(all.counts,data.frame(i = j[2]))
}

colnames(all.counts) <- samples

save(all.counts,file='all.counts.mg.Rdata')
View(all.counts)

merged.counts <- all.counts[-(55358:55362),]
View(merged.counts)
tail(merged.counts)

ensembl.all = read.table('gencode.vM28.annotation.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');

save(ensembl.code, file = "ensembl.code.Rdata")

rm(ensembl.all)
rm(ensembl.gene.names)
rm(ensembl.id.names)

load("ensembl.code.Rdata")

merged.counts = merge(merged.counts, ensembl.code, by="row.names", all.x=F)
rownames(merged.counts) <- make.names(merged.counts$gene, unique=TRUE)
merged.counts <- merged.counts[,-c(1,66,67)]


rm(ensembl.code)
View(merged.counts)

save(merged.counts, file="merged.counts.mg.Rdata")
load("merged.counts.mg.Rdata")

sample.conditions = factor(c(rep("MG_PolyIC_E12_F",4), 
                             rep("MG_PolyIC_E12_M",4), 
                             rep("MG_PolyIC_E14_F",4), 
                             rep("MG_PolyIC_E14_M",4),
                             rep("MG_PolyIC_P5_F",4), 
                             rep("MG_PolyIC_P5_M",4), 
                             rep("MG_PolyIC_P90_F",4), 
                             rep("MG_PolyIC_P90_M",4),
                             rep("MG_Saline_E12_F",4), 
                             rep("MG_Saline_E12_M",4), 
                             rep("MG_Saline_E14_F",4), 
                             rep("MG_Saline_E14_M",4),
                             rep("MG_Saline_P5_F",4), 
                             rep("MG_Saline_P5_M",4), 
                             rep("MG_Saline_P90_F",4), 
                             rep("MG_Saline_P90_M",4)))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions);
dds = DESeq(deseq.counts.table)
vsd = varianceStabilizingTransformation(dds)
plotPCA(vsd)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic() 

view(pcaData)
pcaData2<-separate(pcaData,
                          col=condition,
                          into = c("Tissue","Treatment","Timepoint","Sex"),
                          sep = "_",
                          remove = TRUE)
view(pcaData2)

tiff("mg_pca_sex_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaData2, aes(PC1, PC2, color=Sex, shape=Timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()

tiff("mg_pca_tx_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaData2, aes(PC1, PC2, color=Treatment, shape=Timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","dodgerblue4"))+
  theme_classic() 
dev.off()

ggplot(pcaData2, aes(PC1, PC2, color=Treatment, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","blue"))+
  labs(title = "Microglia") +
  theme_classic() 




## next cluster within each time point!

#E12
View(merged.counts)
mergedE12<-merged.counts[,c(1:8,33:40)]
sample.conditions = factor(c(rep("MG_PolyIC_E12_F",4), 
                             rep("MG_PolyIC_E12_M",4),
                             rep("MG_Saline_E12_F",4), 
                             rep("MG_Saline_E12_M",4)))

deseq.counts.table.e12 = DESeqDataSetFromMatrix(mergedE12, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table.e12)$condition<-factor(colData(deseq.counts.table.e12)$sample.conditions);
dds = DESeq(deseq.counts.table.e12)
vsd = vst(dds)
plotPCA(vsd)
pcaDataE12 <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaDataE12, "percentVar"))

pcaDataE12<-separate(pcaDataE12,
                   col=condition,
                   into = c("Tissue","Treatment","Timepoint","Sex"),
                   sep = "_",
                   remove = TRUE)

tiff("mg_pca_e12_tx_sex.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataE12, aes(PC1, PC2, color=Treatment, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","dodgerblue4"))+
  theme_classic() 
dev.off()

tiff("mg_pca_e12_sex_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataE12, aes(PC1, PC2, color=Sex, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()

#E14
View(merged.counts)
mergedE14<-merged.counts[,c(9:16,41:48)]
sample.conditions = factor(c(rep("MG_PolyIC_E14_F",4), 
                             rep("MG_PolyIC_E14_M",4),
                             rep("MG_Saline_E14_F",4), 
                             rep("MG_Saline_E14_M",4)))

deseq.counts.table.e14 = DESeqDataSetFromMatrix(mergedE14, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table.e14)$condition<-factor(colData(deseq.counts.table.e14)$sample.conditions);
dds = DESeq(deseq.counts.table.e14)
vsd = vst(dds)
plotPCA(vsd)
pcaDataE14 <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaDataE14, "percentVar"))

pcaDataE14<-separate(pcaDataE14,
                     col=condition,
                     into = c("Tissue","Treatment","Timepoint","Sex"),
                     sep = "_",
                     remove = TRUE)

tiff("mg_pca_e14_tx_sex.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataE14, aes(PC1, PC2, color=Treatment, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","dodgerblue4"))+
  theme_classic() 
dev.off()

tiff("mg_pca_e14_sex_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataE14, aes(PC1, PC2, color=Sex, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()

#P5
View(merged.counts)
mergedP5<-merged.counts[,c(17:24,49:56)]
sample.conditions = factor(c(rep("MG_PolyIC_P5_F",4), 
                             rep("MG_PolyIC_P5_M",4),
                             rep("MG_Saline_P5_F",4), 
                             rep("MG_Saline_P5_M",4)))

deseq.counts.table.p5 = DESeqDataSetFromMatrix(mergedP5, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table.p5)$condition<-factor(colData(deseq.counts.table.p5)$sample.conditions);
dds = DESeq(deseq.counts.table.p5)
vsd = vst(dds)
plotPCA(vsd)
pcaDataP5 <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaDataP5, "percentVar"))

pcaDataP5<-separate(pcaDataP5,
                     col=condition,
                     into = c("Tissue","Treatment","Timepoint","Sex"),
                     sep = "_",
                     remove = TRUE)

tiff("mg_pca_p5_tx_sex.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataP5, aes(PC1, PC2, color=Treatment, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","dodgerblue4"))+
  theme_classic() 
dev.off()

tiff("mg_pca_p5_sex_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataP5, aes(PC1, PC2, color=Sex, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()


#P90
View(merged.counts)
mergedP90<-merged.counts[,c(25:32,57:64)]
sample.conditions = factor(c(rep("MG_PolyIC_P90_F",4), 
                             rep("MG_PolyIC_P90_M",4),
                             rep("MG_Saline_P90_F",4), 
                             rep("MG_Saline_P90_M",4)))

deseq.counts.table.p90 = DESeqDataSetFromMatrix(mergedP90, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table.p90)$condition<-factor(colData(deseq.counts.table.p90)$sample.conditions);
dds = DESeq(deseq.counts.table.p90)
vsd = vst(dds)
plotPCA(vsd)
pcaDataP90 <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaDataP90, "percentVar"))

pcaDataP90<-separate(pcaDataP90,
                    col=condition,
                    into = c("Tissue","Treatment","Timepoint","Sex"),
                    sep = "_",
                    remove = TRUE)

tiff("mg_pca_p90_tx_sex.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataP90, aes(PC1, PC2, color=Treatment, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkorange","dodgerblue4"))+
  theme_classic() 
dev.off()

tiff("mg_pca_p90_sex_timepoint.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataP90, aes(PC1, PC2, color=Sex, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()


## SALINE ONLY - all time points

view(merged.counts)
mergedSaline<-merged.counts[,c(33:64)]
View(mergedSaline)
sample.conditions = factor(c(rep("MG_Saline_E12_F",4), 
                             rep("MG_Saline_E12_M",4),
                             rep("MG_Saline_E14_F",4),
                             rep("MG_Saline_E14_M",4),
                             rep("MG_Saline_P5_F",4), 
                             rep("MG_Saline_P5_M",4),
                             rep("MG_Saline_P90_F",4),
                             rep("MG_Saline_P90_M",4)))

deseq.counts.table.saline = DESeqDataSetFromMatrix(mergedSaline, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table.saline)$condition<-factor(colData(deseq.counts.table.saline)$sample.conditions);
dds = DESeq(deseq.counts.table.saline)
vsd = vst(dds)
plotPCA(vsd)
pcaDataSaline <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaDataSaline, "percentVar"))

pcaDataSaline<-separate(pcaDataSaline,
                     col=condition,
                     into = c("Tissue","Treatment","Timepoint","Sex"),
                     sep = "_",
                     remove = TRUE)

tiff("mg_pca_saline_timepoint_sex.tiff", units="in", width=4, height=3, res=300)
ggplot(pcaDataSaline, aes(PC1, PC2, color=Timepoint, shape=Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("darkgreen","darkseagreen","skyblue","dodgerblue4"))+
  theme_classic() 
dev.off()







####Differential expression analysis #####

#### MG_PolyIC_E12 F vs M####

load('merged.counts.mg.Rdata')

View(merged.counts)
merged.counts.small = merged.counts[,c(1:8)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E12_F",4), rep("MG_PolyIC_E12_M",4)), levels=c("MG_PolyIC_E12_F","MG_PolyIC_E12_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = as.data.frame(results(mm.atac))

MG_PolyIC_E12_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_PolyIC_E12_FvsM.lattice

write.csv(MG_PolyIC_E12_FvsM.lattice, file = "MG_PolyIC_E12_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG_PolyIC_E14 F vs M####
merged.counts
merged.counts.small = merged.counts[,9:16]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E14_F",4), rep("MG_PolyIC_E14_M",4)), levels=c("MG_PolyIC_E14_F","MG_PolyIC_E14_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_PolyIC_E14_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_PolyIC_E14_FvsM.lattice

write.csv(MG_PolyIC_E14_FvsM.lattice, file = "MG_PolyIC_E14_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG_PolyIC_P5 F vs M####
View(merged.counts)
merged.counts.small = merged.counts[,c(17:24)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P5_F",4), rep("MG_PolyIC_P5_M",4)), levels=c("MG_PolyIC_P5_F","MG_PolyIC_P5_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_PolyIC_P5_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_PolyIC_P5_FvsM.lattice

write.csv(MG_PolyIC_P5_FvsM.lattice, file = "MG_PolyIC_P5_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG_PolyIC_P90 F vs M####
View(merged.counts)
merged.counts.small = merged.counts[,c(25:32)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P90_F",4), rep("MG_PolyIC_P90_M",4)), levels=c("MG_PolyIC_P90_F","MG_PolyIC_P90_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_PolyIC_P90_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')

MG_PolyIC_P90_FvsM.lattice

write.csv(MG_PolyIC_P90_FvsM.lattice, file = "MG_PolyIC_P90_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG_Saline_E12 F vs M####


View(merged.counts)
merged.counts.small = merged.counts[,c(33:40)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_E12_F",4), rep("MG_Saline_E12_M",4)), levels=c("MG_Saline_E12_F","MG_Saline_E12_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = as.data.frame(results(mm.atac))

MG_Saline_E12_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_Saline_E12_FvsM.lattice

write.csv(MG_Saline_E12_FvsM.lattice, file = "MG_Saline_E12_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG_PolyIC_E14 F vs M####
merged.counts
merged.counts.small = merged.counts[,41:48]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_E14_F",4), rep("MG_Saline_E14_M",4)), levels=c("MG_Saline_E14_F","MG_Saline_E14_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_Saline_E14_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_Saline_E14_FvsM.lattice

write.csv(MG_Saline_E14_FvsM.lattice, file = "MG_Saline_E14_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG_PolyIC_P5 F vs M####
View(merged.counts)
merged.counts.small = merged.counts[,c(49:56)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_P5_F",4), rep("MG_Saline_P5_M",4)), levels=c("MG_Saline_P5_F","MG_Saline_P5_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_Saline_P5_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')
MG_Saline_P5_FvsM.lattice

write.csv(MG_Saline_P5_FvsM.lattice, file = "MG_Saline_P5_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG_PolyIC_P90 F vs M####
View(merged.counts)
merged.counts.small = merged.counts[,c(57:64)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_P90_F",4), rep("MG_Saline_P90_M",4)), levels=c("MG_Saline_P90_F","MG_Saline_P90_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
normalize.counts.5xFADv5xFADSykKO <- counts(mm.atac,normalized = TRUE)
res.mm.atac = results(mm.atac)

MG_Saline_P90_FvsM.lattice =
  categorize.deseq.df(res.mm.atac,
                      thresh = 0.1, log2fold = 0.0, treat = 'M')

MG_Saline_P90_FvsM.lattice

write.csv(MG_Saline_P90_FvsM.lattice, file = "MG_Saline_P90_MvsF.lattice_041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG PolyIC vs Saline E12 F ####
View(merged.counts)
merged.counts.small = merged.counts[,c(1:4,33:36)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E12_F",4), rep("MG_Saline_E12_F",4)), levels=c("MG_PolyIC_E12_F","MG_Saline_E12_F"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E12_F.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_E12_F.lattice

write.csv(MG_PolyIC_vs_Saline_E12_F.lattice, file = "MG_PolyIC_vs_Saline_E12_F.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)




#### MG PolyIC vs Saline E14 F ####

View(merged.counts)
merged.counts.small = merged.counts[,c(9:12,41:44)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E14_F",4), rep("MG_Saline_E14_F",4)), levels=c("MG_PolyIC_E14_F","MG_Saline_E14_F"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E14_F.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_E14_F.lattice

write.csv(MG_PolyIC_vs_Saline_E14_F.lattice, file = "MG_PolyIC_vs_Saline_E14_F.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)

#### MG PolyIC vs Saline P5 F ####

View(merged.counts)
merged.counts.small = merged.counts[,c(17:20,49:52)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P5_F",4), rep("MG_Saline_P5_F",4)), levels=c("MG_PolyIC_P5_F","MG_Saline_P5_F"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P5_F.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_P5_F.lattice

write.csv(MG_PolyIC_vs_Saline_P5_F.lattice, file = "MG_PolyIC_vs_Saline_P5_F.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)

#### MG PolyIC vs Saline P90 F ####

View(merged.counts)
merged.counts.small = merged.counts[,c(25:28,57:60)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P90_F",4), rep("MG_Saline_P90_F",4)), levels=c("MG_PolyIC_P90_F","MG_Saline_P90_F"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P90_F.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_P90_F.lattice

write.csv(MG_PolyIC_vs_Saline_P90_F.lattice, file = "MG_PolyIC_vs_Saline_P90_F.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG PolyIC vs Saline E12 M ####
View(merged.counts)
merged.counts.small = merged.counts[,c(5:8,37:40)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E12_M",4), rep("MG_Saline_E12_M",4)), levels=c("MG_PolyIC_E12_M","MG_Saline_E12_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E12_M.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_E12_M.lattice

write.csv(MG_PolyIC_vs_Saline_E12_M.lattice, file = "MG_PolyIC_vs_Saline_E12_M.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)




#### MG PolyIC vs Saline E14 M ####

View(merged.counts)
merged.counts.small = merged.counts[,c(13:16,45:48)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E14_M",4), rep("MG_Saline_E14_M",4)), levels=c("MG_PolyIC_E14_M","MG_Saline_E14_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E14_M.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_E14_M.lattice

write.csv(MG_PolyIC_vs_Saline_E14_M.lattice, file = "MG_PolyIC_vs_Saline_E14_M.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG PolyIC vs Saline P5 F ####

View(merged.counts)
merged.counts.small = merged.counts[,c(21:24,53:56)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P5_M",4), rep("MG_Saline_P5_M",4)), levels=c("MG_PolyIC_P5_M","MG_Saline_P5_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P5_M.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_P5_M.lattice

write.csv(MG_PolyIC_vs_Saline_P5_M.lattice, file = "MG_PolyIC_vs_Saline_P5_M.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG PolyIC vs Saline P90 F ####

View(merged.counts)
merged.counts.small = merged.counts[,c(29:32,61:64)]
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P90_M",4), rep("MG_Saline_P90_M",4)), levels=c("MG_PolyIC_P90_M","MG_Saline_P90_M"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P90_M.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

MG_PolyIC_vs_Saline_P90_M.lattice

write.csv(MG_PolyIC_vs_Saline_P90_M.lattice, file = "MG_PolyIC_vs_Saline_P90_M.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)




#### MG PolyIC vs Saline E12 - both sexes ####

View(merged.counts)
merged.counts.small = merged.counts[,c(1:8,33:40)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E12",8), rep("MG_Saline_E12",8)), levels=c("MG_PolyIC_E12","MG_Saline_E12"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E12_MnF.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

head(MG_PolyIC_vs_Saline_E12_MnF.lattice)
view(MG_PolyIC_vs_Saline_E12_MnF.lattice)

write.csv(MG_PolyIC_vs_Saline_E12_MnF.lattice, file = "MG_PolyIC_vs_Saline_E12_MnF.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG PolyIC vs Saline E14 - both sexes ####

View(merged.counts)
merged.counts.small = merged.counts[,c(9:16,41:48)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_E14",8), rep("MG_Saline_E14",8)), levels=c("MG_PolyIC_E14","MG_Saline_E14"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_E14_MnF.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

head(MG_PolyIC_vs_Saline_E14_MnF.lattice)
view(MG_PolyIC_vs_Saline_E14_MnF.lattice)

write.csv(MG_PolyIC_vs_Saline_E14_MnF.lattice, file = "MG_PolyIC_vs_Saline_E14_MnF.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)



#### MG PolyIC vs Saline P5 - both sexes ####

View(merged.counts)
merged.counts.small = merged.counts[,c(17:24,49:56)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P5",8), rep("MG_Saline_P5",8)), levels=c("MG_PolyIC_P5","MG_Saline_P5"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P5_MnF.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

head(MG_PolyIC_vs_Saline_P5_MnF.lattice)
view(MG_PolyIC_vs_Saline_P5_MnF.lattice)

write.csv(MG_PolyIC_vs_Saline_P5_MnF.lattice, file = "MG_PolyIC_vs_Saline_P5_MnF.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)


#### MG PolyIC vs Saline P90 - both sexes ####

View(merged.counts)
merged.counts.small = merged.counts[,c(25:32,57:64)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_PolyIC_P90",8), rep("MG_Saline_P90",8)), levels=c("MG_PolyIC_P90","MG_Saline_P90"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_PolyIC_vs_Saline_P90_MnF.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'PolyIC')

head(MG_PolyIC_vs_Saline_P90_MnF.lattice)
view(MG_PolyIC_vs_Saline_P90_MnF.lattice)

write.csv(MG_PolyIC_vs_Saline_P90_MnF.lattice, file = "MG_PolyIC_vs_Saline_P90_MnF.lattice.041522.csv", sep = "\t",
          row.names = TRUE, col.names = TRUE)







## MG - saline - E12 to E14 comparison

View(merged.counts)
merged.counts.small = merged.counts[,c(33:40,41:48)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_E12",8), rep("MG_Saline_E14",8)), levels=c("MG_Saline_E12","MG_Saline_E14"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_Saline_MnF_E12_v_E14.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'E14')

head(MG_Saline_MnF_E12_v_E14.lattice)
view(MG_Saline_MnF_E12_v_E14.lattice)

write.csv(MG_Saline_MnF_E12_v_E14.lattice, file = "MG_Saline_MnF_E12_v_E14.lattice.042022.csv",row.names = TRUE)



## MG - saline - E14 to P5 comparison

View(merged.counts)
merged.counts.small = merged.counts[,c(41:48,49:56)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_E14",8), rep("MG_Saline_P5",8)), levels=c("MG_Saline_E14","MG_Saline_P5"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_Saline_MnF_E14_v_P5.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'P5')

head(MG_Saline_MnF_E14_v_P5.lattice)
view(MG_Saline_MnF_E14_v_P5.lattice)

write.csv(MG_Saline_MnF_E14_v_P5.lattice, file = "MG_Saline_MnF_E14_v_P5.lattice.042022.csv",row.names = TRUE)



## MG - saline - P5 to P90 comparison

View(merged.counts)
merged.counts.small = merged.counts[,c(49:56,57:64)]
colnames(merged.counts.small)<-gsub("_F.*","",colnames(merged.counts.small))
colnames(merged.counts.small)<-gsub("_M.*","",colnames(merged.counts.small))
View(merged.counts.small)

sample.conditions = factor(c(rep("MG_Saline_P5",8), rep("MG_Saline_P90",8)), levels=c("MG_Saline_P5","MG_Saline_P90"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)
mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)
x = as.data.frame(res.mm.atac)

MG_Saline_MnF_P5_v_P90.lattice =
  categorize.deseq.df(x,
                      thresh = 0.1, log2fold = 0.0, treat = 'P90')

head(MG_Saline_MnF_P5_v_P90.lattice)
view(MG_Saline_MnF_P5_v_P90.lattice)

write.csv(MG_Saline_MnF_P5_v_P90.lattice, file = "MG_Saline_MnF_P5_v_P90.lattice.042022.csv",row.names = TRUE)












## Save all the lattices
save(MG_PolyIC_E12_FvsM.lattice, file = "MG_PolyIC_E12_FvsM.lattice.Rdata")
save(MG_PolyIC_E14_FvsM.lattice, file = "MG_PolyIC_E14_FvsM.lattice.Rdata")
save(MG_Saline_E12_FvsM.lattice, file = "MG_Saline_E12_FvsM.lattice.Rdata")
save(MG_Saline_E14_FvsM.lattice, file = "MG_Saline_E14_FvsM.lattice.Rdata")
save(MG_PolyIC_P5_FvsM.lattice, file = "MG_PolyIC_P5_FvsM.lattice.Rdata")
save(MG_PolyIC_P90_FvsM.lattice, file = "MG_PolyIC_P90_FvsM.lattice.Rdata")
save(MG_Saline_P5_FvsM.lattice, file = "MG_Saline_P5_FvsM.lattice.Rdata")
save(MG_Saline_P90_FvsM.lattice, file = "MG_Saline_P90_FvsM.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E12_F.lattice, file = "MG_PolyIC_vs_Saline_E12_F.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E14_F.lattice, file = "MG_PolyIC_vs_Saline_E14_F.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E12_M.lattice, file = "MG_PolyIC_vs_Saline_E12_M.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E14_M.lattice, file = "MG_PolyIC_vs_Saline_E14_M.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P5_F.lattice, file = "MG_PolyIC_vs_Saline_P5_F.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P90_F.lattice, file = "MG_PolyIC_vs_Saline_P90_F.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P5_M.lattice, file = "MG_PolyIC_vs_Saline_P5_M.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P90_M.lattice, file = "MG_PolyIC_vs_Saline_P90_M.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E12_MnF.lattice, file = "MG_PolyIC_vs_Saline_E12_MnF.lattice.Rdata")
save(MG_PolyIC_vs_Saline_E14_MnF.lattice, file = "MG_PolyIC_vs_Saline_E14_MnF.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P5_MnF.lattice, file = "MG_PolyIC_vs_Saline_P5_MnF.lattice.Rdata")
save(MG_PolyIC_vs_Saline_P90_MnF.lattice, file = "MG_PolyIC_vs_Saline_P90_MnF.lattice.Rdata")
save(MG_Saline_MnF_E12_v_E14.lattice, file = "MG_Saline_MnF_E12_v_E14.lattice.Rdata")
save(MG_Saline_MnF_E14_v_P5.lattice, file = "MG_Saline_MnF_E14_v_P5.lattice.Rdata")
save(MG_Saline_MnF_P5_v_P90.lattice, file = "MG_Saline_MnF_P5_v_P90.lattice.Rdata")




#### # Differentially Expressed Genes ####

# load the datasets
load("MG_PolyIC_E12_FvsM.lattice.Rdata")
load("MG_PolyIC_E14_FvsM.lattice.Rdata")
load("MG_Saline_E12_FvsM.lattice.Rdata")
load("MG_Saline_E14_FvsM.lattice.Rdata")
load("MG_PolyIC_P5_FvsM.lattice.Rdata")
load("MG_PolyIC_P90_FvsM.lattice.Rdata")
load("MG_Saline_P5_FvsM.lattice.Rdata")
load("MG_Saline_P90_FvsM.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E12_F.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E14_F.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E12_M.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E14_M.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P5_F.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P90_F.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P5_M.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P90_M.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E12_MnF.lattice.Rdata")
load("MG_PolyIC_vs_Saline_E14_MnF.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P5_MnF.lattice.Rdata")
load("MG_PolyIC_vs_Saline_P90_MnF.lattice.Rdata")
load("MG_Saline_MnF_E12_v_E14.lattice.Rdata")
load("MG_Saline_MnF_E14_v_P5.lattice.Rdata")
load("MG_Saline_MnF_P5_v_P90.lattice.Rdata")

load("merged.counts.mg.Rdata")


#show the number of differentially expressed genes
print(table(MG_PolyIC_E12_FvsM.lattice$response))
print(table(MG_PolyIC_E14_FvsM.lattice$response))
print(table(MG_Saline_E12_FvsM.lattice$response))
print(table(MG_Saline_E14_FvsM.lattice$response))
print(table(MG_PolyIC_P5_FvsM.lattice$response))
print(table(MG_PolyIC_P90_FvsM.lattice$response))
print(table(MG_Saline_P5_FvsM.lattice$response))
print(table(MG_Saline_P90_FvsM.lattice$response))
print(table(MG_PolyIC_vs_Saline_E12_F.lattice$response))
print(table(MG_PolyIC_vs_Saline_E14_F.lattice$response))
print(table(MG_PolyIC_vs_Saline_E12_M.lattice$response))
print(table(MG_PolyIC_vs_Saline_E14_M.lattice$response))
print(table(MG_PolyIC_vs_Saline_P5_F.lattice$response))
print(table(MG_PolyIC_vs_Saline_P90_F.lattice$response))
print(table(MG_PolyIC_vs_Saline_P5_M.lattice$response))
print(table(MG_PolyIC_vs_Saline_P90_M.lattice$response))
print(table(MG_PolyIC_vs_Saline_E12_MnF.lattice$response))
print(table(MG_PolyIC_vs_Saline_E14_MnF.lattice$response))
print(table(MG_PolyIC_vs_Saline_P5_MnF.lattice$response))
print(table(MG_PolyIC_vs_Saline_P90_MnF.lattice$response))
print(table(MG_Saline_MnF_E12_v_E14.lattice$response))
print(table(MG_Saline_MnF_E14_v_P5.lattice$response))
print(table(MG_Saline_MnF_P5_v_P90.lattice$response))










#### DEG PLOTS ####
#we will make heatmaps, MA plots, and volcano plots now


# Microglia E12 polyI:C v. saline - combined males and females

view(merged.counts)
polyicF = merged.counts[1:4]
polyicM = merged.counts[5:8]
salineF = merged.counts[33:36]
salineM = merged.counts[37:40]
c = cbind(salineM,salineF,polyicM,polyicF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E12",8), rep("MG_PolyIC_E12",8)), levels=c("MG_Saline_E12","MG_PolyIC_E12"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12","MG_PolyIC_E12"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("")

tiff("volcano_mg_e12_polyic_v_saline.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,300),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




# Microglia E14 polyI:C v. saline - combined males and females

view(merged.counts)
polyicF = merged.counts[9:12]
polyicM = merged.counts[13:16]
salineF = merged.counts[41:44]
salineM = merged.counts[45:48]
c = cbind(salineM,salineF,polyicM,polyicF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E14",8), rep("MG_PolyIC_E14",8)), levels=c("MG_Saline_E14","MG_PolyIC_E14"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E14","MG_PolyIC_E14"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("")

tiff("volcano_mg_e14_polyic_v_saline.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,50),
                xlim = c(-6,6),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()


# Microglia P5 polyI:C v. saline - combined males and females

view(merged.counts)
polyicF = merged.counts[17:20]
polyicM = merged.counts[21:24]
salineF = merged.counts[49:52]
salineM = merged.counts[53:56]
c = cbind(salineM,salineF,polyicM,polyicF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_P5",8), rep("MG_PolyIC_P5",8)), levels=c("MG_Saline_P5","MG_PolyIC_P5"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P5","MG_PolyIC_P5"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Safb","Plekho2","Slc25a25")

tiff("volcano_mg_p5_polyic_v_saline.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,50),
                xlim = c(-6,6),
                labSize = 3,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()


# Microglia P90 polyI:C v. saline - combined males and females

view(merged.counts)
polyicF = merged.counts[25:28]
polyicM = merged.counts[29:32]
salineF = merged.counts[57:60]
salineM = merged.counts[61:64]
c = cbind(salineM,salineF,polyicM,polyicF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_P90",8), rep("MG_PolyIC_P90",8)), levels=c("MG_Saline_P90","MG_PolyIC_P90"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P90","MG_PolyIC_P90"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("")

tiff("volcano_mg_p90_polyic_v_saline.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,50),
                xlim = c(-6,6),
                labSize = 3,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()














# Microglia E12 saline - male v female

view(merged.counts)
salineF = merged.counts[33:36]
salineM = merged.counts[37:40]
c = cbind(salineF,salineM)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E12_M",4), rep("MG_Saline_E12_F",4)), levels=c("MG_Saline_E12_M","MG_Saline_E12_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12_M","MG_Saline_E12_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Kdm5d","Uty","Ddx3y","Eif2s3y","Xist","Eif2s3x","Gm29650")

tiff("volcano_mg_e12_saline_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,80),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()





# Microglia E14 saline - male v female

view(merged.counts)
salineF = merged.counts[41:44]
salineM = merged.counts[45:48]
c = cbind(salineF,salineM)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E14_M",4), rep("MG_Saline_E14_F",4)), levels=c("MG_Saline_E14_M","MG_Saline_E14_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E14_M","MG_Saline_E14_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Gm29650")

tiff("volcano_mg_e14_saline_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,80),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




# Microglia P5 saline - male v female

view(merged.counts)
salineF = merged.counts[49:52]
salineM = merged.counts[53:56]
c = cbind(salineF,salineM)
view(c)

sample.conditions = factor(c(rep("MG_Saline_P5_M",4), rep("MG_Saline_P5_F",4)), levels=c("MG_Saline_P5_M","MG_Saline_P5_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P5_M","MG_Saline_P5_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Xist","Ddx3y","Kdm5d","Uty","Eif2s3y")

tiff("volcano_mg_p5_saline_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,250),
                xlim = c(-9,9),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




# Microglia P90 saline - male v female

view(merged.counts)
salineF = merged.counts[57:60]
salineM = merged.counts[61:64]
c = cbind(salineF,salineM)
view(c)

sample.conditions = factor(c(rep("MG_Saline_P90_M",4), rep("MG_Saline_P90_F",4)), levels=c("MG_Saline_P90_M","MG_Saline_P90_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P90_M","MG_Saline_P90_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Ddx3y","Eif2s3y","Kdm5d","Xist","Uty")

tiff("volcano_mg_p90_saline_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,250),
                xlim = c(-9,9),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()










# Microglia E12 polyic - male v female

view(merged.counts)
polyicF = merged.counts[1:4]
polyicM = merged.counts[5:8]
c = cbind(polyicF,polyicM)
view(c)

sample.conditions = factor(c(rep("MG_PolyIC_E12_M",4), rep("MG_PolyIC_E12_F",4)), levels=c("MG_PolyIC_E12_M","MG_PolyIC_E12_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_PolyIC_E12_M","MG_PolyIC_E12_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Uty","Ddx3y","Xist","Kdm5d","Eif2s3y")

tiff("volcano_mg_e12_polyic_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,80),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




# Microglia E14 polyic - male v female

view(merged.counts)
polyicF = merged.counts[9:12]
polyicM = merged.counts[13:16]
c = cbind(polyicF,polyicM)
view(c)

sample.conditions = factor(c(rep("MG_PolyIC_E14_M",4), rep("MG_PolyIC_E14_F",4)), levels=c("MG_PolyIC_E14_M","MG_PolyIC_E14_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_PolyIC_E14_M","MG_PolyIC_E14_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Uty","Ddx3y","Xist","Kdm5d","Eif2s3y","Gm29650")

tiff("volcano_mg_e14_polyic_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,80),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()





# Microglia P5 polyic - male v female

view(merged.counts)
polyicF = merged.counts[17:20]
polyicM = merged.counts[21:24]
c = cbind(polyicF,polyicM)
view(c)

sample.conditions = factor(c(rep("MG_PolyIC_P5_M",4), rep("MG_PolyIC_P5_F",4)), levels=c("MG_PolyIC_P5_M","MG_PolyIC_P5_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_PolyIC_P5_M","MG_PolyIC_P5_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Uty","Ddx3y","Xist","Kdm5d","Eif2s3y")

tiff("volcano_mg_p5_polyic_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,200),
                xlim = c(-7,7),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()





# Microglia P90 polyic - male v female

view(merged.counts)
polyicF = merged.counts[25:28]
polyicM = merged.counts[29:32]
c = cbind(polyicF,polyicM)
view(c)

sample.conditions = factor(c(rep("MG_PolyIC_P90_M",4), rep("MG_PolyIC_P90_F",4)), levels=c("MG_PolyIC_P90_M","MG_PolyIC_P90_F"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_PolyIC_P90_M","MG_PolyIC_P90_F"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Uty","Ddx3y","Xist","Kdm5d","Eif2s3y")

tiff("volcano_mg_p90_polyic_m_v_f.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-1,250),
                xlim = c(-8,8),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




# Microglia saline - combined males and females - E12 v E14

view(merged.counts)
E12salineF = merged.counts[33:36]
E12salineM = merged.counts[37:40]
E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
c = cbind(E12salineM,E12salineF,E14salineM,E14salineF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E12",8), rep("MG_Saline_E14",8)), levels=c("MG_Saline_E12","MG_Saline_E14"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12","MG_Saline_E14"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Ptgds","Kel","Mmp9","Matn1","Col11a2","Acan","Cytl1","Meltf","Ibsp","Tectb","Onecut3")

tiff("volcano_mg_saline_e12_v_e14.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-5,100),
                xlim = c(-6,6),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34', 'thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()





# Microglia saline - combined males and females - E14 v P5

view(merged.counts)
E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
c = cbind(E14salineM,E14salineF,P5salineM,P5salineF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_E14",8), rep("MG_Saline_P5",8)), levels=c("MG_Saline_E14","MG_Saline_P5"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E14","MG_Saline_P5"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Col2a1","Tgfa","Cbr2","Mlxipl","P3h2","Col9a1","Mest","Acan","Rnase2a","Retnla","Cd209g","Slc47a1","Crybb1","P2ry12","Sall1","Egr1","Mef2a","Jun","Dcx","Sox14","Neurog1")

tiff("volcano_mg_saline_e14_v_p5.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-5,300),
                xlim = c(-10,10),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34','thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()





# Microglia saline - combined males and females - P5 v P90

view(merged.counts)
P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
P90salineF = merged.counts[57:60]
P90salineM = merged.counts[61:64]
c = cbind(P5salineM,P5salineF,P90salineM,P90salineF)
view(c)

sample.conditions = factor(c(rep("MG_Saline_P5",8), rep("MG_Saline_P90",8)), levels=c("MG_Saline_P5","MG_Saline_P90"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P5","MG_Saline_P90"));
dds = DESeq(deseq.counts.table)
coefList <- resultsNames(dds)
resLFC.vol <- lfcShrink(dds, coef = 2, type="apeglm")

geneLabels <- c("Sting1","Selplg","Tmem119","Lpl","Spp1","Dab2","Mcm5","Csf1","P2ry12","Psat1","Cx3cr1","Csf1r","Sall1","Pmepa1")

tiff("volcano_mg_saline_p5_v_p90.tiff", units="in", width=8, height=6, res=300)
EnhancedVolcano(resLFC.vol,
                ylim = c(-5,320),
                xlim = c(-10,10),
                labSize = 4,
                lab= rownames(resLFC.vol),
                drawConnectors = TRUE,
                endsConnectors = 'last',
                max.overlaps = 100,
                lengthConnectors = unit(.01,"npc"),
                selectLab = geneLabels,
                x='log2FoldChange',
                y='padj',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                col=c('gray', 'gray34','thistle', 'mediumorchid4'),
                subtitle = NULL,
                title = NULL)
dev.off()




## HEATMAPS

# MG Saline - time point comparisons

# E12 v E14

E12salineF = merged.counts[33:36]
E12salineM = merged.counts[37:40]
E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
c = cbind(E12salineM,E12salineF,E14salineM,E14salineF)
sample.conditions = factor(c(rep("MG_Saline_E12",8), rep("MG_Saline_E14",8)), levels=c("MG_Saline_E12","MG_Saline_E14"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12","MG_Saline_E14"));
dds = DESeq(deseq.counts.table)

lattice = MG_Saline_MnF_E12_v_E14.lattice
lattice = na.omit(lattice)
rld_HH = rlogTransformation(dds)

# top 10 UP and DOWN genes

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:10]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:10]

a = assay(rld_HH)[c(inc[1:10],dec[1:10]),]
a = a - rowMeans(a)

dev.off()
tiff("heatmap_mg_saline_MnF_E12_v_E14_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","darkseagreen","darkgreen"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()

dev.off()
tiff("heatmap_mg_saline_MnF_E12_v_E14_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","thistle","mediumorchid4"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()

dev.off()
tiff("heatmap_mg_saline_MnF_E12_v_E14_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("mediumorchid4","thistle","whitesmoke", 
                            "whitesmoke","darkseagreen","darkgreen"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()




# E14 v P5

E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
c = cbind(E14salineM,E14salineF,P5salineM,P5salineF)
sample.conditions = factor(c(rep("MG_Saline_E14",8), rep("MG_Saline_P5",8)), levels=c("MG_Saline_E14","MG_Saline_P5"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E14","MG_Saline_P5"));
dds = DESeq(deseq.counts.table)

lattice = MG_Saline_MnF_E14_v_P5.lattice
lattice = na.omit(lattice)
rld_HH = rlogTransformation(dds)

# top 10 UP and DOWN genes

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:10]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:10]

a = assay(rld_HH)[c(inc[1:10],dec[1:10]),]
a = a - rowMeans(a)

dev.off()
tiff("heatmap_mg_saline_MnF_E14_v_P5_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","darkseagreen","darkgreen"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()

dev.off()
tiff("heatmap_mg_saline_MnF_E14_v_P5_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","thistle","mediumorchid4"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()



# P5 v P90

P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
P90salineF = merged.counts[57:60]
P90salineM = merged.counts[61:64]
c = cbind(P5salineM,P5salineF,P90salineM,P90salineF)
sample.conditions = factor(c(rep("MG_Saline_P5",8), rep("MG_Saline_P90",8)), levels=c("MG_Saline_P5","MG_Saline_P90"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P5","MG_Saline_P90"));
dds = DESeq(deseq.counts.table)

lattice = MG_Saline_MnF_P5_v_P90.lattice
lattice = na.omit(lattice)
rld_HH = rlogTransformation(dds)

# top 10 UP and DOWN genes

y = lattice[lattice$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:10]

y = lattice[lattice$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:10]

a = assay(rld_HH)[c(inc[1:10],dec[1:10]),]
a = a - rowMeans(a)

dev.off()
tiff("heatmap_mg_saline_MnF_P5_v_P90_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","darkseagreen","darkgreen"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()

dev.off()
tiff("heatmap_mg_saline_MnF_P5_v_P90_top20degs.tiff", units="in", width=4, height=3.5, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","thistle","mediumorchid4"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()




# key microglial genes expressed over time

E12salineF = merged.counts[33:36]
E12salineM = merged.counts[37:40]
E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
P90salineF = merged.counts[57:60]
P90salineM = merged.counts[61:64]
c = cbind(E12salineM,E12salineF,E14salineM,E14salineF,P5salineM,P5salineF,P90salineM,P90salineF)
view(c)


sample.conditions = factor(c(rep("MG_Saline_E12",8),rep("MG_Saline_E14",8),rep("MG_Saline_P5",8), rep("MG_Saline_P90",8)), levels=c("MG_Saline_E12","MG_Saline_E14","MG_Saline_P5","MG_Saline_P90"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12","MG_Saline_E14","MG_Saline_P5","MG_Saline_P90"));
dds = DESeq(deseq.counts.table)
rld_HH = vst(dds)


a<-c("Dab2","Mcm5","Pf4","F13a1","Cdk1","Scd2","Psat1","Csf1","Crybb1","Sall1","Egr1","Fcrls","Cd14","Mef2a","Irf8","Cx3cr1","Jun","Fos","Sting1","Mafb","P2ry12","Tmem119")

a = assay(rld_HH)[c(a),]
a = a - rowMeans(a)

dev.off()
tiff("heatmap_mg_saline_MnF_timecourse.tiff", units="in", width=5, height=4, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","darkseagreen","darkgreen"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()

dev.off()
tiff("heatmap_mg_saline_MnF_timecourse.tiff", units="in", width=5, height=4, res=300)
pheatmap(a, 
         colorRampPalette(c("dodgerblue4","skyblue3","whitesmoke", 
                            "whitesmoke","thistle","mediumorchid4"))(100), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         show_rownames=T, 
         show_colnames = FALSE)
dev.off()













### Pathway analyses ####

# MG - saline only - males and females combined




# E12 v E14

E12salineF = merged.counts[33:36]
E12salineM = merged.counts[37:40]
E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
c = cbind(E12salineM,E12salineF,E14salineM,E14salineF)
sample.conditions = factor(c(rep("MG_Saline_E12",8), rep("MG_Saline_E14",8)), levels=c("MG_Saline_E12","MG_Saline_E14"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E12","MG_Saline_E14"));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = MG_Saline_MnF_E12_v_E14.lattice
lattice = na.omit(lattice)

res <- results(dds)
resOrdered <- res[order(res$padj),] 
resSig <- subset(resOrdered, padj < 0.1) 
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)
upGenes.list <- rownames(resSigUp)
downGenes.list <- rownames(resSigDown)

gostres <- gost(query = list("Upregulated Genes" = upGenes.list, 
                             "Downregulated Genes" = downGenes.list),
                organism = "mmusculus", 
                ordered_query = TRUE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = TRUE,
                user_threshold = 0.05, correction_method = "gSCS", 
                domain_scope = "annotated",
                sources = c("GO", "KEGG", "REAC", "TF", "CORUM", "MIRNA", "WP"),
                as_short_link = FALSE)

gostres_link <- gost(list("Upregulated Genes" = upGenes.list,
                          "Downregulated Genes" = downGenes.list), 
                     as_short_link = TRUE)
print(gostres_link)




# E14 v P5

E14salineF = merged.counts[41:44]
E14salineM = merged.counts[45:48]
P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
c = cbind(E14salineM,E14salineF,P5salineM,P5salineF)
sample.conditions = factor(c(rep("MG_Saline_E14",8), rep("MG_Saline_P5",8)), levels=c("MG_Saline_E14","MG_Saline_P5"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_E14","MG_Saline_P5"));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = MG_Saline_MnF_E14_v_P5.lattice
lattice = na.omit(lattice)

res <- results(dds)
resOrdered <- res[order(res$padj),] 
resSig <- subset(resOrdered, padj < 0.1) 
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)
upGenes.list <- rownames(resSigUp)
downGenes.list <- rownames(resSigDown)

gostres <- gost(query = list("Upregulated Genes" = upGenes.list, 
                             "Downregulated Genes" = downGenes.list),
                organism = "mmusculus", 
                ordered_query = TRUE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = TRUE,
                user_threshold = 0.05, correction_method = "gSCS", 
                domain_scope = "annotated",
                sources = c("GO", "KEGG", "REAC", "TF", "CORUM", "MIRNA", "WP"),
                as_short_link = FALSE)

gostres_link <- gost(list("Upregulated Genes" = upGenes.list,
                          "Downregulated Genes" = downGenes.list), 
                     as_short_link = TRUE)
print(gostres_link)





# P5 v P90

P5salineF = merged.counts[49:52]
P5salineM = merged.counts[53:56]
P90salineF = merged.counts[57:60]
P90salineM = merged.counts[61:64]
c = cbind(P5salineM,P5salineF,P90salineM,P90salineF)
sample.conditions = factor(c(rep("MG_Saline_P5",8), rep("MG_Saline_P90",8)), levels=c("MG_Saline_P5","MG_Saline_P90"))
deseq.counts.table = DESeqDataSetFromMatrix(c, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$sample.conditions<-factor(colData(deseq.counts.table)$sample.conditions, levels=c("MG_Saline_P5","MG_Saline_P90"));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

lattice = MG_Saline_MnF_P5_v_P90.lattice
lattice = na.omit(lattice)

res <- results(dds)
resOrdered <- res[order(res$padj),] 
resSig <- subset(resOrdered, padj < 0.1) 
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigDown <- subset(resSig, log2FoldChange < 0)
upGenes.list <- rownames(resSigUp)
downGenes.list <- rownames(resSigDown)

gostres <- gost(query = list("Upregulated Genes" = upGenes.list, 
                             "Downregulated Genes" = downGenes.list),
                organism = "mmusculus", 
                ordered_query = TRUE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = TRUE,
                user_threshold = 0.05, correction_method = "gSCS", 
                domain_scope = "annotated",
                sources = c("GO", "KEGG", "REAC", "TF", "CORUM", "MIRNA", "WP"),
                as_short_link = FALSE)

gostres_link <- gost(list("Upregulated Genes" = upGenes.list,
                          "Downregulated Genes" = downGenes.list), 
                     as_short_link = TRUE)
print(gostres_link)






##### reactome term plots! #####



## e12 v e14

df=gProfiler_MG_Saline_MnF_E12_v_E14
reac.terms.e12<-df[df$source == 'REAC', ]


reac.terms<-as.data.frame(reac.terms.e12$term_name)
names(reac.terms)[1]<-'Term'
reac.terms$TermSize<-reac.terms.e12$term_size

reac.terms$GenesUp<-reac.terms.e12$`query_size__Upregulated Genes`
reac.terms$GenesDown<-reac.terms.e12$`query_size__Downregulated Genes`

reac.terms$IntersectionSizeUp<-reac.terms.e12$`intersection_size__Upregulated Genes`
reac.terms$IntersectionSizeUp[reac.terms$IntersectionSizeUp<0] <-- NA
reac.terms$IntersectionSizeDown<-reac.terms.e12$`intersection_size__Downregulated Genes`
reac.terms$IntersectionSizeDown[reac.terms$IntersectionSizeDown<0] <-- NA

reac.terms$GeneRatioUp<-reac.terms$IntersectionSizeUp/reac.terms$TermSize
reac.terms$GeneRatioDown<-reac.terms$IntersectionSizeDown/reac.terms$TermSize

reac.terms<-reac.terms %>%
  mutate(GeneRatio = 
           GeneRatioUp %>% 
           is.na %>%
           ifelse(GeneRatioDown, GeneRatioUp) )

reac.terms$padjUp<-reac.terms.e12$`adjusted_p_value__Upregulated Genes`
reac.terms$padjDown<-reac.terms.e12$`adjusted_p_value__Downregulated Genes`
reac.terms$log10adjpUp <-- log10(reac.terms$padjUp)
reac.terms$log10adjpDown <-- log10(reac.terms$padjDown)
reac.terms$log10adjP<-- reac.terms$log10adjpUp + reac.terms$log10adjpDown
reac.terms$log10adjpUp[reac.terms$log10adjpUp==0] <-- NA
reac.terms$log10adjpDown[reac.terms$log10adjpDown==0] <-- NA

reac.terms$log10adjP<--0-reac.terms$log10adjP #swap direction of p value for correct graphing below

view(reac.terms)

reac.terms.handpicked.e12<-reac.terms[c(1,3,13,31,39,49,70,99,109,151),]
View(reac.terms.handpicked)

reac.terms.handpicked.e12 <- reac.terms.handpicked.e12 %>% 
  mutate(Color = ifelse(log10adjpDown > 0, "darkgreen",
                        ifelse(log10adjpUp > 0, "darkseagreen")))
reac.terms.handpicked.e12$Color<-replace_na(reac.terms.handpicked.e12$Color, "darkseagreen")



dev.off()
tiff("reac.mg.saline.e12.v.e14.dot.tiff", units="in", width=6.5, height=6.5, res=300)
ggplot(reac.terms.handpicked.e12, 
       aes(x = GeneRatio, 
           y =fct_reorder(Term, GeneRatio))) + 
  geom_point(aes(size = abs(log10adjP), color = Color))+
  scale_color_identity() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(size = 10, colour = "black")) +
  ylab(NULL)  
dev.off()



rm(reac.terms.handpicked.e12)
reac.terms.handpicked.e12<-reac.terms[c(1,3,13,31,39,49,70,99,109,151),]

reac.terms.handpicked.e12 <- reac.terms.handpicked.e12 %>% 
  mutate(Timepoint = ifelse(log10adjpDown > 0, "E12",
                        ifelse(log10adjpUp > 0, "E14")))
reac.terms.handpicked.e12$Timepoint<-replace_na(reac.terms.handpicked.e12$Timepoint, "E14")

reac.terms.handpicked.e12$Term[reac.terms.handpicked.e12$Term == "Cell Cycle"] <- "Cell cycle"
reac.terms.handpicked.e12$Term[reac.terms.handpicked.e12$Term == "DNA Replication"] <- "DNA replication"
reac.terms.handpicked.e12$Term[reac.terms.handpicked.e12$Term == "Mitotic Spindle Checkpoint"] <- "Mitotic spindle checkpoint"
reac.terms.handpicked.e12$Term[reac.terms.handpicked.e12$Term == "DNA Repair"] <- "DNA repair"


View(reac.terms.handpicked.e12)


dev.off()
tiff("reac.mg.saline.e12.v.e14.bar.tiff", units="in", width=9, height=4, res=300)
ggplot(reac.terms.handpicked.e12,
       aes(x = log10adjP,
           y = fct_rev(fct_reorder(Term, Timepoint)), 
           fill= Timepoint))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkgreen","darkseagreen"))+
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(size = 15, colour = "black")) +
  ylab(NULL)
dev.off()





# e14 v p5

df=gProfiler_MG_Saline_MnF_E14_v_P5
reac.terms.e14<-df[df$source == 'REAC', ]


reac.terms<-as.data.frame(reac.terms.e14$term_name)
names(reac.terms)[1]<-'Term'
reac.terms$TermSize<-reac.terms.e14$term_size

reac.terms$GenesUp<-reac.terms.e14$`query_size__Upregulated Genes`
reac.terms$GenesDown<-reac.terms.e14$`query_size__Downregulated Genes`

reac.terms$IntersectionSizeUp<-reac.terms.e14$`intersection_size__Upregulated Genes`
reac.terms$IntersectionSizeUp[reac.terms$IntersectionSizeUp<0] <-- NA
reac.terms$IntersectionSizeDown<-reac.terms.e14$`intersection_size__Downregulated Genes`
reac.terms$IntersectionSizeDown[reac.terms$IntersectionSizeDown<0] <-- NA

reac.terms$GeneRatioUp<-reac.terms$IntersectionSizeUp/reac.terms$TermSize
reac.terms$GeneRatioDown<-reac.terms$IntersectionSizeDown/reac.terms$TermSize

reac.terms<-reac.terms %>%
  mutate(GeneRatio = 
           GeneRatioUp %>% 
           is.na %>%
           ifelse(GeneRatioDown, GeneRatioUp) )

reac.terms$padjUp<-reac.terms.e14$`adjusted_p_value__Upregulated Genes`
reac.terms$padjDown<-reac.terms.e14$`adjusted_p_value__Downregulated Genes`
reac.terms$log10adjpUp <-- log10(reac.terms$padjUp)
reac.terms$log10adjpDown <-- log10(reac.terms$padjDown)
reac.terms$log10adjP<-- reac.terms$log10adjpUp + reac.terms$log10adjpDown
reac.terms$log10adjpUp[reac.terms$log10adjpUp==0] <-- NA
reac.terms$log10adjpDown[reac.terms$log10adjpDown==0] <-- NA

reac.terms$log10adjP<--0-reac.terms$log10adjP #swap direction of p value for correct graphing below

view(reac.terms)

reac.terms.handpicked.e14<-reac.terms[c(),]
View(reac.terms.handpicked.e14)

reac.terms.handpicked.e14 <- reac.terms.handpicked.e14 %>% 
  mutate(Timepoint = ifelse(log10adjpDown > 0, "E14",
                            ifelse(log10adjpUp > 0, "P5")))
reac.terms.handpicked.e14$Timepoint<-replace_na(reac.terms.handpicked.e14$Timepoint, "P5")

reac.terms.handpicked.e14$Term[reac.terms.handpicked.e14$Term == "Potassium Channels"] <- "Potassium channels"
reac.terms.handpicked.e14$Term[reac.terms.handpicked.e14$Term == "Death Receptor Signlaing"] <- "Death receptor signaling"
reac.terms.handpicked.e14$Term[reac.terms.handpicked.e14$Term == "Toll-like Receptor Cascades"] <- "Toll-like receptor cascades"
reac.terms.handpicked.e14$Term[reac.terms.handpicked.e14$Term == "Collagen chain trimerization"] <- " Collagen chain trimerization"

View(reac.terms.handpicked.e14)


dev.off()
tiff("reac.mg.saline.e14.v.p5.bar.tiff", units="in", width=9, height=4, res=300)
ggplot(reac.terms.handpicked.e14,
       aes(x = log10adjP,
           y = fct_rev(fct_reorder(Term, Timepoint)), 
           fill= Timepoint))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkseagreen","skyblue"))+
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(size = 15, colour = "black")) +
  ylab(NULL)
dev.off()






# p5 to p90

df=gProfiler_MG_Saline_MnF_P5_v_P90
reac.terms.p5<-df[df$source == 'REAC', ]


reac.terms<-as.data.frame(reac.terms.p5$term_name)
names(reac.terms)[1]<-'Term'
reac.terms$TermSize<-reac.terms.p5$term_size

reac.terms$GenesUp<-reac.terms.p5$`query_size__Upregulated Genes`
reac.terms$GenesDown<-reac.terms.p5$`query_size__Downregulated Genes`

reac.terms$IntersectionSizeUp<-reac.terms.p5$`intersection_size__Upregulated Genes`
reac.terms$IntersectionSizeUp[reac.terms$IntersectionSizeUp<0] <-- NA
reac.terms$IntersectionSizeDown<-reac.terms.p5$`intersection_size__Downregulated Genes`
reac.terms$IntersectionSizeDown[reac.terms$IntersectionSizeDown<0] <-- NA

reac.terms$GeneRatioUp<-reac.terms$IntersectionSizeUp/reac.terms$TermSize
reac.terms$GeneRatioDown<-reac.terms$IntersectionSizeDown/reac.terms$TermSize

reac.terms<-reac.terms %>%
  mutate(GeneRatio = 
           GeneRatioUp %>% 
           is.na %>%
           ifelse(GeneRatioDown, GeneRatioUp) )

reac.terms$padjUp<-reac.terms.p5$`adjusted_p_value__Upregulated Genes`
reac.terms$padjDown<-reac.terms.p5$`adjusted_p_value__Downregulated Genes`
reac.terms$log10adjpUp <-- log10(reac.terms$padjUp)
reac.terms$log10adjpDown <-- log10(reac.terms$padjDown)
reac.terms$log10adjP<-- reac.terms$log10adjpUp + reac.terms$log10adjpDown
reac.terms$log10adjpUp[reac.terms$log10adjpUp==0] <-- NA
reac.terms$log10adjpDown[reac.terms$log10adjpDown==0] <-- NA

reac.terms$log10adjP<--0-reac.terms$log10adjP #swap direction of p value for correct graphing below

view(reac.terms)

reac.terms.handpicked.p5<-reac.terms[c(1,7,18,38,57,70,77,108,115,140),]
View(reac.terms.handpicked.p5)

reac.terms.handpicked.p5 <- reac.terms.handpicked.p5 %>% 
  mutate(Timepoint = ifelse(log10adjpDown > 0, "P5",
                            ifelse(log10adjpUp > 0, "P90")))
reac.terms.handpicked.p5$Timepoint<-replace_na(reac.terms.handpicked.p5$Timepoint, "P90")

reac.terms.handpicked.p5$Term[reac.terms.handpicked.p5$Term == "Axon guidance"] <- " Axon guidance"
reac.terms.handpicked.p5$Term[reac.terms.handpicked.p5$Term == "Cell Cycle"] <- " Cell cycle"
reac.terms.handpicked.p5$Term[reac.terms.handpicked.p5$Term == "Cytokine Signaling in Immune system"] <- "Cytokine signaling in immune system"
reac.terms.handpicked.p5$Term[reac.terms.handpicked.p5$Term == "RNA Polymerase II Transcription"] <- "RNA polymerase II transcription"


View(reac.terms.handpicked.p5)


dev.off()
tiff("reac.mg.saline.p5.v.p90.bar.tiff", units="in", width=9, height=4, res=300)
ggplot(reac.terms.handpicked.p5,
       aes(x = log10adjP,
           y = fct_rev(fct_reorder(Term, Timepoint)), 
           fill= Timepoint))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("skyblue","dodgerblue4"))+
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(size = 15, colour = "black")) +
  ylab(NULL)
dev.off()




## combine all together into a single plot

reac.summary <- dplyr::bind_rows(list(reac.terms.handpicked.e12, reac.terms.handpicked.e14,reac.terms.handpicked.p5), .id = 'source')
view(reac.summary)

reac.summary$source[reac.summary$source == "1"] <- "E12 v. E14"
reac.summary$source[reac.summary$source == "2"] <- "P5 v. E14"
reac.summary$source[reac.summary$source == "3"] <- "P90 v. P5"
names(reac.summary)[1]<-"Category"

terms<-reac.summary$Term
view(terms)

reac.summary$Term <- factor(reac.summary$Term, levels = c("Cell cycle","Metabolism of RNA","DNA replication","Mitotic spindle checkpoint", "DNA repair", "Collagen chain trimerization", "ECM proteoglycans", "Degradation of the extracellular matrix", "Collagen degradation","Integrin cell surface interactions","Protein-protein interactions at synapses", "Extracellular matrix organization", " Collagen chain trimerization","Potassium channels","Axon guidance","Metabolism of lipids","Neutrophil degranulation","Toll-like receptor cascades","Death Receptor Signalling","TNF signaling"," Cell cycle","Translation","Synthesis of DNA","rRNA processing"," Axon guidance","RNA polymerase II transcription","DAP12 signaling","Cytokine signaling in immune system","Signaling by SCF-KIT", "Class I MHC mediated antigen processing & presentation"))



dev.off()
tiff("reac.mg.saline.reac.summary.bar.tiff", units="in", width=12, height=7, res=300)
ggplot(reac.summary,
       aes(x = log10adjP,
           y = fct_rev(fct_reorder(Term,Category)), 
           fill= Timepoint))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("darkgreen","darkseagreen","skyblue","dodgerblue4"))+
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(size = 15, colour = "black")) +
  ylab(NULL)
dev.off()

# bless up that this finally worked :,)



