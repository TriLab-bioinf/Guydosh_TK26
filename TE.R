#Project goal: to run DEseq and calculate TE data for 2-5A and poly I:C treated samples

library(pacman)
p_load("tidyverse", "DESeq2", "pheatmap", "ggplot2",'EnhancedVolcano', "gmodels")

#load csv files
Ribo_count <- read.csv("~/Desktop/RNase_L_project/TE_calculations/June2022/analysis/data/Ribo_count.csv")
RNA_count <- read.csv("~/Desktop/RNase_L_project/TE_calculations/June2022/analysis/data/RNA_count.csv")


# preparing structure for DEseq
cdsfoot_table<-Ribo_count[,c(7:19)]
cdsfoot_table<-round(cdsfoot_table) #this rounds up or down the fraction numbers
cdsrna_table<-RNA_count[,c(3:15)]
cdsrna_table<-round(cdsrna_table)

# for the fold changes I decided to run these separately

# Ribo-seq only
#this is for WT 2-5A Riboseq only
countTable<-cdsfoot_table[,c(1,6,4,11,12)]
condition=factor(c("condition1","condition1","condition2","condition2", "condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable))
dds<-DESeqDataSetFromMatrix(countTable,coldata,~treatment)
dds_out<-DESeq(dds)
res<-results(dds_out)
rownames(res)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res),file="control_vs_25A_Riboseq.csv")

#this is for WT Riboseq IC only
countTable_IC<-cdsfoot_table[,c(1,6,5, 10)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_IC))
dds_IC<-DESeqDataSetFromMatrix(countTable_IC,coldata,~treatment)
dds_out_IC<-DESeq(dds_IC)
res_IC<-results(dds_out_IC)
rownames(res_IC)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_IC),file="control_vs_IC_Riboseq.csv")

#this is for KO 2-5A Riboseq only
countTable_KO<-cdsfoot_table[,c(2,7,3,8)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_KO))
dds_KO<-DESeqDataSetFromMatrix(countTable_KO,coldata,~treatment)
dds_out_KO<-DESeq(dds_KO)
res_KO<-results(dds_out_KO)
rownames(res_KO)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_KO),file="KOcontrol_vs_25A_Riboseq.csv")

#this is for KO Riboseq IC only
countTable_KO_IC<-cdsfoot_table[,c(2,7,9,13)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_KO_IC))
dds_KO_IC<-DESeqDataSetFromMatrix(countTable_KO_IC,coldata,~treatment)
dds_out_KO_IC<-DESeq(dds_KO_IC)
res_KO_IC<-results(dds_out_KO_IC)
rownames(res_KO_IC)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_KO_IC),file="KOcontrol_vs_IC_Riboseq.csv")

#RNA-seq only
#this is for WT 2-5A RNA-seq only
countTable_RNA<-cdsrna_table[,c(1,6,4,11,12)]
condition=factor(c("condition1","condition1","condition2","condition2", "condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_RNA))
dds_RNA<-DESeqDataSetFromMatrix(countTable_RNA,coldata,~treatment)
dds_out_RNA<-DESeq(dds_RNA)
res_RNA<-results(dds_out_RNA)
rownames(res_RNA)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_RNA),file="control_vs_25A_RNAseq.csv")

#this is for WT RNA-seq IC only
countTable_IC_RNA<-cdsrna_table[,c(1,6,5, 10)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_IC_RNA))
dds_IC_RNA<-DESeqDataSetFromMatrix(countTable_IC_RNA,coldata,~treatment)
dds_out_IC_RNA<-DESeq(dds_IC_RNA)
res_IC_RNA<-results(dds_out_IC_RNA)
rownames(res_IC_RNA)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_IC_RNA),file="control_vs_IC_RNAseq.csv")

#this is for KO 2-5A RNA-seq only
countTable_KO_RNA<-cdsrna_table[,c(2,7,3,8)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_KO_RNA))
dds_KO_RNA<-DESeqDataSetFromMatrix(countTable_KO_RNA,coldata,~treatment)
dds_out_KO_RNA<-DESeq(dds_KO_RNA)
res_KO_RNA<-results(dds_out_KO_RNA)
rownames(res_KO_RNA)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_KO_RNA),file="KOcontrol_vs_25A_RNAseq.csv")

#this is for KO RNAseq IC only
countTable_KO_IC_RNA<-cdsrna_table[,c(2,7,9,13)]
condition=factor(c("condition1","condition1","condition2","condition2"))
coldata=DataFrame(treatment=condition,row.names=colnames(countTable_KO_IC_RNA))
dds_KO_IC_RNA<-DESeqDataSetFromMatrix(countTable_KO_IC_RNA,coldata,~treatment)
dds_out_KO_IC_RNA<-DESeq(dds_KO_IC_RNA)
res_KO_IC_RNA<-results(dds_out_KO_IC_RNA)
rownames(res_KO_IC_RNA)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_KO_IC_RNA),file="KOcontrol_vs_IC_RNAseq.csv")

#TE_calculations
#This is for WT 2-5A treated
TE_countTable<-cbind(cdsfoot_table[,c(1,6,4,11,12)], cdsrna_table[,c(1,6,4,11,12)])
condition=factor(c("condition1","condition1","condition2","condition2", "condition2", "condition1","condition1","condition2","condition2", "condition2"))
type=factor(c("foot","foot","foot","foot","foot","mRNA","mRNA","mRNA","mRNA","mRNA"))
coldata=DataFrame(treatment=condition,exp=type,row.names=colnames(TE_countTable))
dds_TE<-DESeqDataSetFromMatrix(TE_countTable,coldata,~treatment + exp + treatment:exp)
dds_TE$treatment <- relevel(dds_TE$treatment,ref="condition1")
dds_TE$exp <- relevel(dds_TE$exp,ref="mRNA")
dds_out_TE<-DESeq(dds_TE)
res_TE<-results(dds_out_TE)
rownames(res_TE)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_TE),file="TE_WT_controlvs25A.csv")

#vulcano plot
setEPS()
postscript( "TE_WT_control_25A.eps", width=6.5, height=8)
EnhancedVolcano(res_TE,
                lab = rownames(res_TE),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(-0.5, 7),
                title = 'Control vs 2-5A',
                titleLabSize = 20,
                col=c('grey', 'grey', 'grey', 'red3'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                shape = c(1, 1, 1, 1),
                labSize = 3,
                pointSize = 3,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)

dev.off() 

#This is for WT poly I:C treated
TE_IC_countTable<-cbind(cdsfoot_table[,c(1,6,5,10)], cdsrna_table[,c(1,6,5,10)])
condition=factor(c("condition1","condition1","condition2","condition2", "condition1","condition1","condition2", "condition2"))
type=factor(c("foot","foot","foot","foot","mRNA","mRNA","mRNA","mRNA"))
coldata=DataFrame(treatment=condition,exp=type,row.names=colnames(TE_IC_countTable))
dds_TE_IC<-DESeqDataSetFromMatrix(TE_IC_countTable,coldata,~treatment + exp + treatment:exp)
dds_TE_IC$treatment <- relevel(dds_TE_IC$treatment,ref="condition1")
dds_TE_IC$exp <- relevel(dds_TE_IC$exp,ref="mRNA")
dds_out_TE_IC<-DESeq(dds_TE_IC)
res_TE_IC<-results(dds_out_TE_IC)
rownames(res_TE_IC)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_TE_IC),file="TE_WT_controlvspolyIC.csv")

#vulcano plot
setEPS()
postscript( "TE_WT_control_polyIC.eps", width=6.5, height=8)
EnhancedVolcano(res_TE_IC,
                lab = rownames(res_TE_IC),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(-0.5, 7),
                title = 'Control vs poly I:C',
                titleLabSize = 20,
                col=c('grey', 'grey', 'grey', 'red3'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                shape = c(1, 1, 1, 1),
                labSize = 3,
                pointSize = 3,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)

dev.off() 

#This is for KO 2-5A treated
TE_KO_countTable<-cbind(cdsfoot_table[,c(2,7,3,8)], cdsrna_table[,c(2,7,3,8)])
condition=factor(c("condition1","condition1","condition2", "condition2", "condition1","condition1","condition2", "condition2"))
type=factor(c("foot","foot","foot","foot","mRNA","mRNA","mRNA","mRNA"))
coldata=DataFrame(treatment=condition,exp=type,row.names=colnames(TE_KO_countTable))
dds_TE_KO<-DESeqDataSetFromMatrix(TE_KO_countTable,coldata,~treatment + exp + treatment:exp)
dds_TE_KO$treatment <- relevel(dds_TE_KO$treatment,ref="condition1")
dds_TE_KO$exp <- relevel(dds_TE_KO$exp,ref="mRNA")
dds_out_TE_KO<-DESeq(dds_TE_KO)
res_TE_KO<-results(dds_out_TE_KO)
rownames(res_TE_KO)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_TE_KO),file="TE_KO_controlvs25A.csv")

#vulcano plot
setEPS()
postscript( "TE_KO_control_25A.eps", width=6.5, height=8)
EnhancedVolcano(res_TE_KO,
                lab = rownames(res_TE_KO),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(-0.5, 7),
                title = 'Control vs 2-5A RNase L KO',
                titleLabSize = 20,
                col=c('grey', 'grey', 'grey', 'red3'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                shape = c(1, 1, 1, 1),
                labSize = 3,
                pointSize = 3,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)

dev.off() 

#This is for KO poly I:C treated
TE_IC_KO_countTable<-cbind(cdsfoot_table[,c(2,7,9,13)], cdsrna_table[,c(2,7,9,13)])
condition=factor(c("condition1","condition1","condition2","condition2", "condition1","condition1","condition2", "condition2"))
type=factor(c("foot","foot","foot","foot","mRNA","mRNA","mRNA","mRNA"))
coldata=DataFrame(treatment=condition,exp=type,row.names=colnames(TE_IC_KO_countTable))
dds_TE_IC_KO<-DESeqDataSetFromMatrix(TE_IC_KO_countTable,coldata,~treatment + exp + treatment:exp)
dds_TE_IC_KO$treatment <- relevel(dds_TE_IC_KO$treatment,ref="condition1")
dds_TE_IC_KO$exp <- relevel(dds_TE_IC_KO$exp,ref="mRNA")
dds_out_TE_IC_KO<-DESeq(dds_TE_IC_KO)
res_TE_IC_KO<-results(dds_out_TE_IC_KO)
rownames(res_TE_IC_KO)<-RNA_count$Gene_symbol
write.csv(as.data.frame(res_TE_IC_KO),file="TE_WT_controlvspolyIC.csv")

#vulcano plot
setEPS()
postscript( "TE_KO_control_polyIC.eps", width=6.5, height=8)
EnhancedVolcano(res_TE_IC_KO,
                lab = rownames(res_TE_IC_KO),
                x = 'log2FoldChange',
                y = 'pvalue',
                ylim = c(-0.5, 7),
                title = 'Control vs poly I:C RNase L KO',
                titleLabSize = 20,
                col=c('grey', 'grey', 'grey', 'red3'),
                pCutoff = 0.05,
                FCcutoff = 1.5,
                shape = c(1, 1, 1, 1),
                labSize = 3,
                pointSize = 3,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE
)

dev.off() 

#GO analysis
library(clusterProfiler)
library(enrichplot)

organism = "org.Hs.eg.db"
library(org.Hs.eg.db)

#preparing data for 2-5A
df = read.csv("TE_WT_controlvs25A.csv", header=TRUE)
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$X
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

require(DOSE)
setEPS()
postscript( "GO_TE_WT_controlvs25A.eps")
dotplot(gse, showCategory=10, split=".sign", font.size=8) + facet_grid(.~.sign)
dev.off()

#preparing data for poly I:C
df_IC = read.csv("TE_WT_controlvspolyIC.csv", header=TRUE)
original_gene_list_IC <- df_IC$log2FoldChange
names(original_gene_list_IC) <- df_IC$X
gene_list_IC<-na.omit(original_gene_list_IC)
gene_list_IC = sort(gene_list_IC, decreasing = TRUE)
gse_IC <- gseGO(geneList=gene_list_IC, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

setEPS()
postscript( "GO_TE_WT_controlvspolyIC.eps")
dotplot(gse_IC, showCategory=10, split=".sign", font.size=8) + facet_grid(.~.sign)
dev.off()


# creating Venn-diagrams for TEs

p_load(VennDiagram)
library(VennDiagram)

## Data preparation
pval_threshold <- 0.05
WT_25A.degs <- row.names(res_TE[which(res_TE$pvalue <= pval_threshold), ])
WT_IC.degs <- row.names(res_TE_IC[which(res_TE_IC$pvalue <= pval_threshold), ])
KO_25A.degs <- row.names(res_TE_KO[which(res_TE_KO$pvalue <= pval_threshold), ])
KO_IC.degs <- row.names(res_TE_IC_KO[which(res_TE_IC_KO$pvalue <= pval_threshold), ])


venn.plot2 <- draw.quad.venn(length(KO_25A.degs),length(WT_25A.degs),
                                length(KO_IC.degs),length(WT_IC.degs),
                                # Calculate the intersection of the two sets
                                length( intersect(KO_25A.degs, WT_25A.degs) ),
                                length( intersect(KO_25A.degs, KO_IC.degs) ),
                                length( intersect(KO_25A.degs, WT_IC.degs) ),
                                length( intersect(WT_25A.degs, KO_IC.degs) ),
                                length( intersect(WT_25A.degs, WT_IC.degs) ),
                                length( intersect(KO_IC.degs, WT_IC.degs) ),
                                length( intersect( intersect(KO_25A.degs, WT_25A.degs),intersect(WT_25A.degs, KO_IC.degs))),
                                length( intersect(intersect(KO_25A.degs, WT_25A.degs),intersect(WT_25A.degs, WT_IC.degs) )),
                                length( intersect(intersect(KO_25A.degs, KO_IC.degs),intersect(KO_IC.degs, WT_IC.degs) )),
                                length( intersect(intersect(WT_25A.degs, KO_IC.degs),intersect(KO_IC.degs, WT_IC.degs) )),
                                length( intersect(intersect(KO_25A.degs, WT_25A.degs),intersect(KO_IC.degs, WT_IC.degs)) ),
                                category = c("KO_25A", "WT_25A", "KO_IC","WT_IC"),
                                fill=c(" light blue", "light green", "pink", "light yellow"))

# Actually plot the plot
setEPS()
postscript( "all_Venn.eps")
grid.draw(venn.plot2)
dev.off() 
