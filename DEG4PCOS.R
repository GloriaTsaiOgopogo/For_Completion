################## 1 处理表达值（芯片、rpkm等，小数数值）#############

#if(interactive())
if(!require("BiocManager")) install.packages('BiocManager')

######## Reading GSE series matrix:
# Example: GSE137684_series_matrix.txt

cat("Select an data matrix:\n")
path <- file.choose() 
rawdata <- readLines(path, n = -1)
# unlink(path) # tidy up
rawdata <- rawdata[which(rawdata!="")] 
rawdata <- rawdata[!grepl('^!', rawdata)] 
rawmatrix <- read.table(text = rawdata, header = T, sep = "\t")

######## the first column must be the probe names or gene annotations:
tags <- rawmatrix[,1]  

######## Define sample groups:
# colnum number of 2 groups:
# control group: 3,7,9,10
# experimental group: 2,4,5,6,8,11,12,13

tmp <- readline("Define control group: ")
group1 <- as.integer(array(strsplit(tmp,",")[[1]]))
group1_exp <- rawmatrix[,group1]

cat("group1 sample names: \n")
print(colnames(group1_exp))

tmp <- readline("Define experiment group: ")
group2 <- as.integer(array(strsplit(tmp,",")[[1]]))
group2_exp <- rawmatrix[,group2]

cat("group2 sample names: \n")
print(colnames(group2_exp))


if(!require("limma")) BiocManager::install('limma', force = TRUE)
if(!require("stringr")) BiocManager::install('stringr', force = TRUE)

library(limma)
library(stringr)

exp <- cbind(group1_exp,group2_exp)
rownames(exp) <- tags
#boxplot(exp)		

group_list=c(rep('ctr',ncol(group1_exp)),rep("exper",ncol(group2_exp)))
group_list <- factor(group_list,levels = c("ctr","exper"))	
boxplot(exp, outline=FALSE, notch=T, col=group_list, las=2)

## check if normalization needed
exp <- log(exp + 1)
#exp <- log2(exp + 1)

## batch normalize
exp <- normalizeBetweenArrays(exp)
design <- model.matrix(~0+factor(group_list))
colnames(design)=c('ctr','pcos')
#design <- model.matrix(~factor(group_list))

fit <- lmFit(exp,design)
cont.matrix=makeContrasts('ctr-pcos',levels = design)
fit <- contrasts.fit(fit,cont.matrix)
fit <- eBayes(fit)
#options(digits = 6)
deg <- topTable(fit,adjust='BH',number = Inf) %>% na.omit()
head(deg) 

deg[which(deg$logFC >= 1 & deg$adj.P.Val < 0.01),'direction'] <- 'up'
deg[which(deg$logFC <= -1 & deg$adj.P.Val  < 0.01),'direction'] <- 'down'
#deg[which(abs(deg$logFC) <= 1 | deg$adj.P.Val  >= 0.05),'direction'] <- 'none'
summary(factor(deg$direction))
deg <- deg[order(deg$adj.P.Val), ]

write.table(deg, file.choose(), row.names = T, col.names = T, sep = "\t")

	
################## 1 处理counts（整数数值）#############

if(!require("DESeq2")) BiocManager::install('DESeq2', force = TRUE)
library(DESeq2)
library(stringr)

######## Reading GSE series matrix:
# Example: GSE155489_gc_pcos_counts.csv

cat("Select an data matrix:\n")
path <- file.choose() 
rawmatrix <- read.table(path, header = F, sep = ",")  ### notice the delimiter ("\t" or ",")
tags <- rawmatrix[,1]

######## Define sample groups:

tmp <- readline("Define control group: ")
#2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
group1 <- as.integer(array(strsplit(tmp,",")[[1]]))
group1_exp <- rawmatrix[,group1]

cat("group1 sample names: \n")
print(colnames(group1_exp))

tmp <- readline("Define experiment group: ")
#25,26,27,28,29,30,31,32,33,34,35,36,37,38,39
group2 <- as.integer(array(strsplit(tmp,",")[[1]]))
group2_exp <- rawmatrix[,group2]

group_list <- c(rep('ctr',ncol(group1_exp)),rep('exper',ncol(group2_exp)))
group_list <- factor(group_list,levels = c("ctr","exper"))

counts <- cbind(group1_exp,group2_exp)
row.names(counts) <- tags

sample <- data.frame(row.names = colnames(counts), condition = group_list)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample, design = ~condition)

# filter low counts: 
cutoff <- readline("counts filter set: ")
dds <- dds[rowSums(counts(dds))>= as.integer(cutoff), ] 

# normalize:
exp <- vst(dds, blind = FALSE) 
plotPCA(exp, intgroup='condition')
# head(as.data.frame(assay(exp)))

dds <- DESeq(dds)
contrast <- c('condition','ctr','exper')
deg <- results(dds,contrast = contrast) %>% na.omit()
deg <- deg[order(deg$padj),]

write.table(deg, file.choose(), row.names = T, col.names = T, sep = "\t")

deg$change = ifelse(deg$padj < 0.01 & abs(deg$log2FoldChange) >= 1, 
                     ifelse(deg$log2FoldChange> 1 ,'Up','Down'),
                     'Stable')

table(deg$change)



ggplot(
  #设置数据
  as.data.frame(deg), 
  aes(x = log2FoldChange, 
      y = -log10(padj), 
      colour = change)) +
      geom_point(alpha=0.4, size=3.5) +
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
      geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
      geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
# 坐标轴
      labs(x="log2(fold change)",
      y="-log10 (p-value)")+
      theme_bw()+

      # 图例
      theme(plot.title = element_text(hjust = 0.5), 
      legend.position="right", 
      legend.title = element_blank()
)

