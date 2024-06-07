#Load packages

library(DESeq2)
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Conv/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("HC",21),rep("Conv",10)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Conv","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Conv/Result.txt",sep="\t", quote=FALSE,col.names = NA)

###
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Mild/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("HC",21),rep("Mild",26)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Mild","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Mild/Result.txt",sep="\t", quote=FALSE,col.names = NA)

##
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Severe/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("HC",21),rep("Severe",11)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Severe","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC-Severe/Result.txt",sep="\t", quote=FALSE,col.names = NA)

##
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Conv-Mild/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Conv",10),rep("Mild",26)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Mild","Conv"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Conv-Mild/Result.txt",sep="\t", quote=FALSE,col.names = NA)

##
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Conv-Severe/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Conv",10),rep("Severe",11)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Severe","Conv"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Conv-Severe/Result.txt",sep="\t", quote=FALSE,col.names = NA)

##
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Mild",26),rep("Severe",11)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Severe","Mild"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/Result.txt",sep="\t", quote=FALSE,col.names = NA)

##
#### Differential Expression 2

# HC_Conv -vs- Hospitalized-mild
# HC_Conv -vs- Hospitalized-severe


count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("HC",31),rep("Mild",26)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Mild","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Result.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("HC",31),rep("Severe",11)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Severe","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/Result.txt",sep="\t", quote=FALSE,col.names = NA)

############

library(umap)
data=read.delim("/home/anoop/Desktop/COVID_Omics/UMAP/Anova_Signi.txt",header=T,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/UMAP/Anova_Signi_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/UMAP/Anova_Signi_UMAP.txt",row.names = 1,check.names = FALSE)
head(dat)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/UMAP/Anova_Signi_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(HC="#38725e",HC_antibody_positive="#4c8ec1",Hospitalized_mild="#cc689c",Hospitalized_severe="#ccb741"))+
  scale_fill_manual(values=c(HC="#478f76",HC_antibody_positive="#5fb2f2",Hospitalized_mild="#ff83c3",Hospitalized_severe="#e5ce49"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=10),legend.position = c(0.86, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()


library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/UMAP/Anova_Signi_UMAP_label.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+geom_text_repel(data=dat,aes(x=V1,y=V2,label = row.names(dat)),
                                                                                                         size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  scale_color_manual(values=c(HC="#38725e",HC_antibody_positive="#4c8ec1",Hospitalized_mild="#cc689c",Hospitalized_severe="#ccb741"))+
  scale_fill_manual(values=c(HC="#478f76",HC_antibody_positive="#5fb2f2",Hospitalized_mild="#ff83c3",Hospitalized_severe="#e5ce49"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=10),legend.position = c(0.86, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/Master.txt",header=T,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/MasterUMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/MasterUMAP.txt",row.names = 1,check.names = FALSE)
head(dat)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/MasterUMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(HC="#335b76",HC_antibody_positive="#9eb4cc",Hospitalized_mild="#ccac00",Hospitalized_severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",HC_antibody_positive="#C6E2FF",Hospitalized_mild="#ffd700",Hospitalized_severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=10),legend.position = c(0.15, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/Master.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/metadata.tab",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/Master_PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Master_PCA1.txt",row.names = 1,check.names = FALSE)

pdf("/home/anoop/Desktop/COVID_Omics/MasterPCA_label.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = Label),
                                                                                                           size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",HC_antibody_positive="#9eb4cc",Hospitalized_mild="#ccac00",Hospitalized_severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",HC_antibody_positive="#C6E2FF",Hospitalized_mild="#ffd700",Hospitalized_severe="#ffa500"))+
  labs(x="PC1, 47.4% variance",y="PC2, 7.73% variance")+theme(axis.title = element_text(size=10),legend.position = c(0.87, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

ip=read.delim("/home/anoop/Desktop/COVID_Omics/Test/Input.tab",check.names = FALSE,header = T,row.names = 1)
head(ip)
library("ggpubr")
library(reshape2)
library(car)
install.packages("RVAideMemoire")
m=melt(ip)
head(m)
ggplot(data=ip, aes(x=IL8 , group=Group, fill=Group)) +
  geom_density(adjust=1.5, alpha=.4)

hist(ip$IL8)
shapiro.test(IL8 ~ Group, data = ip)

######################## HeatMap

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/Data.txt",header = T,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table((X$carpet),file="/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/Zscore.txt",header = T, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/meta.txt",row.names = 1,header = T)


#col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2a4e6c","#4682b4" ,"#fffcc9","#ffa700", "#ff5a00"))

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

col_fun1 = colorRamp2(c(2, 0, -3), c("#2a4e6c","yellow", "#ff5a00"))
library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Antibody_positive"="#C6E2FF","Covid_mild"="#ffd700","Covid_severe"="#ffa500"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/FC_trans.txt",header = T,check.names = FALSE,row.names = 1)

head(LFC)
ncol(Zscore)
H1=Heatmap(as.matrix(t(Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,row_split = 2,show_column_names = FALSE,
           row_dend_width = unit(2, "cm"),column_title_gp =gpar(fontsize = 0),
           top_annotation  =columnAnnotation(Group = sampleinfo$Group,col=colours,show_legend=TRUE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(30, "cm"),width  = unit(18, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))

col_fun_lfc = colorRamp2(c(-1, -0.05,0, 0.5,1), c("#0000cc","#4c4cff" ,"white", "#b2b200","#7f7f00"))

H3=Heatmap(as.matrix(t(LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(5, "cm"),show_row_names = FALSE,
           row_names_gp =gpar(fontsize = 7),height  = unit(3, "cm"),column_names_gp =gpar(fontsize = 12),na_col = "#e6e6e6") 

tt=H1 + H3 
pdf("/home/anoop/Desktop/COVID_Omics/Krushkal/HeatMap/HeatMap1.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()
?rowAnnotation


pval=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/Kruskal_Wallis_ResultCopy.tab",header = T)
head(pval)
pval$AdjPvalKW =p.adjust(pval$KW,method = "BH")
write.table(pval,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/Kruskal_Wallis_ResultCopyADJ.tab",sep="\t",col.names = NA, quote = FALSE)



##################### HC -vs- Covid >Cluster

data=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_Signi.txt",header=T,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniUMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniUMAP.txt",row.names = 1,check.names = FALSE)
head(dat)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniUMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(HC="#335b76",HC_antibody_positive="#9eb4cc",Hospitalized_mild="#ccac00",Hospitalized_severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",HC_antibody_positive="#C6E2FF",Hospitalized_mild="#ffd700",Hospitalized_severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_Signi.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/metadata.tab",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniPCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniPCA.txt",row.names = 1,check.names = FALSE)

pdf("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniPCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",HC_antibody_positive="#9eb4cc",Hospitalized_mild="#ccac00",Hospitalized_severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",HC_antibody_positive="#C6E2FF",Hospitalized_mild="#ffd700",Hospitalized_severe="#ffa500"))+
  labs(x="PC1, 55.2% variance",y="PC2, 7.1% variance")+theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(1,1,1,1, "cm"),
                                                              legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

rm(list = ls())
packageVersion("limma")
library(limma)
ip=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/KEGG/Input.txt",header = T)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/Krushkal/KEGG/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=nGene,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#191978",high="#7f7fb4",breaks=c(1E-45,1E-10,1E-7,1E-6,1E-5,0.001))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2,6,4,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()
install.packages("distill")
################################# NEW ANALYSIS
library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/PCA/Input.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/PCA/metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/PCA/Master_PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/PCA/Master_PCA.txt",row.names = 1,check.names = FALSE)

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/PCA/MasterPCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",HC_anti_pos="#6f9035",hosp_mild="#ccac00",hosp_severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",HC_anti_pos="#8BB443",hosp_mild="#ffd700",hosp_severe="#ffa500"))+
  labs(x="PC1, 47.6% variance",y="PC2, 7.8% variance")+theme(axis.title = element_text(size=13),axis.text = element_text(size=12),legend.position = c(0.81, 0.91),plot.margin = margin(1,1,1,1, "cm"),
                                                              legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()


pval=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MW_Results.tab",header = T)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/MW_ResultsAdjested.tab",sep="\t",col.names = NA, quote = FALSE)




library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Input.txt",header = T,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table((X$carpet),file="/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Zscore.txt",header = T, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/meta.txt",row.names = 1,header = T)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Antibody_positive"="#C6E2FF","Covid_mild"="#ffd700","Covid_severe"="#ffa500"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/LFC.txt",header = T,check.names = FALSE,row.names = 1)

head(LFC)
ncol(t(Zscore))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,row_split = 2,show_column_names = FALSE,
           row_dend_width = unit(2, "cm"),column_title_gp =gpar(fontsize = 0),
           top_annotation  =columnAnnotation(Group = sampleinfo$Group,col=colours,show_legend=TRUE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(30, "cm"),width  = unit(18, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))

col_fun_lfc = colorRamp2(c(-1, -0.05,0, 0.5,1), c("#0000cc","#4c4cff" ,"white", "#b2b200","#7f7f00"))

H3=Heatmap(as.matrix(t(LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(5, "cm"),show_row_names = FALSE,
           right_annotation = ha,row_names_gp =gpar(fontsize = 7),height  = unit(3, "cm"),column_names_gp =gpar(fontsize = 12),na_col = "#e6e6e6") 


ha = rowAnnotation(foo = anno_mark(at = c(1,2,3,4,7,9,10,11,12,13,15,16,17,18,19,21,22,23,25,26,27,28,29,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,55,56,57,58,60,61,62,63,64,65,66,67,68,69,72,73,74,75,77,78,79,80,81,83,86,87,88,89),
                                   labels_gp = gpar(fontsize=7),lines_gp = gpar(col=NA),link_width=unit(0, "mm"),labels = c("ADA","ADGRG1","ANGPT1","ANGPT2","CASP-8","CCL19","CCL20","CCL23","CCL3","CCL4","CD27","CD28","CD4","CD40","CD40-L","CD70","CD83","CD8A","CSF-1","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL12","CXCL13","CXCL9","DCN","EGF","FASLG","FGF2","Gal-1","Gal-9","GZMA","GZMB","GZMH","HGF","HO-1","ICOSLG","IFN-gamma","IL10","IL12","IL12RB1","IL15","IL18","IL33","IL5","IL6","IL7","IL8","KLRD1","LAG3","LAMP3","LAP_TGF-beta-1","MCP-1","MCP-2","MCP-3","MCP-4","MIC-A_B","MMP12","NCR1","NOS3","PD-L1","PD-L2","PDGF_subunit_B","PGF","PTN","TIE2","TNF","TNFRSF21","TNFSF14","TRAIL","TWEAK","VEGFA")))

CORR=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Correlation.txt",header = T,check.names = FALSE,row.names = 1)
head(CORR)

H4=Heatmap(as.matrix((CORR)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(5, "cm"),show_row_names = FALSE,
           right_annotation = ha,row_names_gp =gpar(fontsize = 7),height  = unit(3, "cm"),column_names_gp =gpar(fontsize = 12),na_col = "#e6e6e6") 

tt=H1 + H3 + H4
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/HeatMap/Heatmap.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()



#############################################
library(psych)
Genes=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Covid.txt",check.names=FALSE,header = T)
Res=corr.test(as.matrix(Genes),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Results.txt",sep="\t",col.names = NA,quote = FALSE)


DAT=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Tmp.txt")
head(DAT)
Re <- acast(DAT,Gene1~Gene2)
head(Re)
write.table(Re,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Significant.txt",sep="\t",col.names = NA,quote = FALSE)
install.packages("ggcorrplot")
library(ggcorrplot)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Corrplot.pdf",width = 10,height = 10)
DD=ggcorrplot(Res$r, p.mat=Res$p,sig.level = 0.1,hc.order = TRUE,insig = "blank",lab_col="black",colors=c("#0000cc","white","#cc0000"),
           tl.cex=5,tl.srt=90,legend.title="Correlation",outline.color="#e6e6e6")
dev.off()
?ggcorrplot()
head(DD$data)
write.table(DD$data,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Ordered_1.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Ordered_1.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))
head(Re)
write.table(Re,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Significant.txt",sep="\t",col.names = NA,quote = FALSE)
library(ComplexHeatmap)
library(circlize)
R=rownames(Re)
(as.data.frame(R))
write.table(as.data.frame(R),file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/RR.txt",sep="\t",col.names = NA,quote = FALSE)
col_corr= colorRamp2(c(-1, -0.5,0, 0.25,0.5,0.75,1), c("#23415a","#3f75a2" ,"white","#e3c156","#db8329","#a91f1f","#401919"))
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Heatmap_1.pdf",width = 12,height = 12)
Heatmap(as.matrix(Re),col=col_corr,cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",height  = unit(19, "cm"),width  = unit(19, "cm"),show_row_names = TRUE,border = TRUE,row_title_gp = gpar(fontsize=0),
        row_gap = unit(0, "mm"),row_names_side = "right",row_names_gp =gpar(fontsize = 5),column_names_gp =gpar(fontsize = 5),na_col = "#e6e6e6")
decorate_annotation("pt", slice = 2,{
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = "#ff5a5a", col = "#ff5a5a"), just = "left")
    grid.text(gp = gpar(fontsize=8),paste("Cytokine-cytokine receptor interaction\nChemokine signaling pathway\nIL-17 signaling pathway\nTNF signaling pathway\nToll-like receptor signaling pathway\nHIF-1 signaling pathway"), x = unit(4, "mm"), just = "left")
  })
dev.off()
########################################### MILD and SEVERE
?grid.text
library(psych)
Genes=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Mild.txt",check.names=FALSE,header = T)
Res=corr.test(as.matrix(Genes),use = "pairwise",method="spearman",adjust="BH")


install.packages("ggcorrplot")
library(ggcorrplot)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Mild.pdf",width = 10,height = 10)
DD=ggcorrplot(Res$r,p.mat = Res$p,sig.level = 0.1,hc.order = TRUE,insig = "blank",lab_col="black",colors=c("#0000cc","white","#cc0000"),
              tl.cex=5,tl.srt=90,legend.title="Correlation",outline.color="#e6e6e6")
dev.off()

write.table(DD$data,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/OrderedMild1.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/OrderedMild1.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))

library(ComplexHeatmap)
library(circlize)
col_corr= colorRamp2(c(-1, -0.05,0, 0.5,1), c("#1919ff","#4c4cff" ,"white", "#ff0000","#cc0000"))
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/HeatmapMild1.pdf",width = 10,height = 10)
Heatmap(as.matrix(Re),col=col_corr,cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width  = unit(17, "cm"),show_row_names = TRUE,border = FALSE,
        row_names_gp =gpar(fontsize = 5),height  = unit(17, "cm"),column_names_gp =gpar(fontsize = 5),na_col = "#e6e6e6")
dev.off()

##############################################

Genes=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Severe.txt",check.names=FALSE,header = T)
Res=corr.test(as.matrix(Genes),use = "pairwise",method="spearman",adjust="BH")


install.packages("ggcorrplot")
library(ggcorrplot)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Severe.pdf",width = 10,height = 10)
DD=ggcorrplot(Res$r,p.mat = Res$p,sig.level = 0.1,hc.order = TRUE,insig = "blank",lab_col="black",colors=c("#0000cc","white","#cc0000"),
              tl.cex=5,tl.srt=90,legend.title="Correlation",outline.color="#e6e6e6")
dev.off()
head(Res$r)
write.table(DD$data,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/OrderedSevere1.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/OrderedSevere1.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))

library(ComplexHeatmap)
library(circlize)
col_corr= colorRamp2(c(-1, -0.05,0, 0.5,1), c("#1919ff","#4c4cff" ,"white", "#ff0000","#cc0000"))
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/HeatmapSevere1.pdf",width = 10,height = 10)
Heatmap(as.matrix(Re),col=col_corr,cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width  = unit(17, "cm"),show_row_names = TRUE,border = TRUE,
        row_names_gp =gpar(fontsize = 5),height  = unit(17, "cm"),column_names_gp =gpar(fontsize = 5),na_col = "#e6e6e6")
dev.off()

#####################################################################Pathway

ip=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/KEGG_HC_Covid/Input.txt",header = T)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/KEGG_HC_Covid/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=nGene,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_colour_gradientn(colours=c('#191978','#a3a3c9'),breaks=c(1E-45,1E-20,1E-5,0.001))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2,6,4,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()

###############################################

install.packages("corrr")
library(corrr)
XX=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Test.txt")
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/Ordered_1.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))

head(Re)
network_plot(Re)

####################################

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Violin_Severe_Mild/Data.txt",check.names = FALSE)


P1=ggplot(Bean,aes(x=Group,y=HGF,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+
  labs(y="NPX",title = "HGF")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=9),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P2=ggplot(Bean,aes(x=Group,y=PTN,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "PTN")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 8, by = 4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P3=ggplot(Bean,aes(x=Group,y=CXCL13,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL13")+
  theme_bw()+scale_y_continuous(breaks = seq(6, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P4=ggplot(Bean,aes(x=Group,y=`MCP-3`,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MCP-3")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 10, by = 3))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P5=ggplot(Bean,aes(x=Group,y=IL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "IL12")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 12, by = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P6=ggplot(Bean,aes(x=Group,y=CCL23,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CCL23")+
  theme_bw()+scale_y_continuous(breaks = seq(8, 14, by = 1.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P7=ggplot(Bean,aes(x=Group,y=TWEAK,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(y="NPX",title = "TWEAK")+
  theme_bw()+scale_y_continuous(breaks = seq(6, 12, by = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=9),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P8=ggplot(Bean,aes(x=Group,y=VEGFA,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "VEGFA")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 12, by = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P9=ggplot(Bean,aes(x=Group,y=CXCL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL12")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 2, by = 0.4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P10=ggplot(Bean,aes(x=Group,y=ANGPT2,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "ANGPT2")+
  theme_bw()+scale_y_continuous(breaks = seq(1, 9, by = 1))+
  theme(axis.title.x = element_blank(),axis.ticks.x=element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
                                                                                                        panel.border = element_blank(),legend.position = "none")

P11=ggplot(Bean,aes(x=Group,y=FASLG,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "FASLG")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = c(1.8, 0.5))+scale_y_continuous(breaks = seq(0, 10, by =0.8))
library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Violin_Severe_Mild/Violin.pdf",height = 5,width =10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,nrow = 2,ncol = 6)+theme(plot.margin = margin(2,0.5,2,0.5, "cm"))
dev.off()

##################################  Circular bubble plot correlation and pathway

library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)


edge=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/BubbleCircular/Edge.txt")
connect=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/BubbleCircular/Connection.txt")
vert=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/BubbleCircular/Vertices.txt")
mygraph <- graph_from_data_frame(edge, vertices=vert )
from = match( connect$from, vert$name)
to = match( connect$to, vert$name)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/Correlation/BubbleCircular/Bubble.pdf")
q=ggraph(mygraph, layout = 'dendrogram', circular = TRUE)
q+geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2,width=0.7,aes(color=group))+scale_edge_colour_manual(values =c("blue","red"))+geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=name), size=2, alpha=1)+geom_node_point(aes(filter = leaf,size=value, alpha=0.2,color=group))+scale_colour_manual(values=c("blue","red"))+theme_void()+theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"cm")) + expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
dev.off()
?geom_conn_bundle



pval=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/Shubha/MW_Results.tab",header = T)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/Shubha/MW_ResultsAdj.tab",sep="\t",col.names = NA, quote = FALSE)


######################################### METABOLOME

library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/Data.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("COVID19_NEG", "COVID19_Convalescent", "Hospitalised_Mild","Hospitalised_Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),
                                                                                                           size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#335b76",COVID19_Convalescent="#9eb4cc",Hospitalised_Mild="#ccac00",Hospitalised_Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#48690E",COVID19_Convalescent="#C6E2FF",Hospitalised_Mild="#ffd700",Hospitalised_Severe="#ffa500"))+
  labs(x="PC1, 40.05% variance",y="PC2, 11.7% variance")+theme(axis.title = element_text(size=10),legend.position = c(0.87, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                                              legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=3.7,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#335b76",COVID19_Convalescent="#9eb4cc",Hospitalised_Mild="#ccac00",Hospitalised_Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#48690E",COVID19_Convalescent="#C6E2FF",Hospitalised_Mild="#ffd700",Hospitalised_Severe="#ffa500"))+
  labs(x="PC1, 40.05% variance",y="PC2, 11.7% variance")+theme(axis.title = element_text(size=10),legend.position = c(0.87, 0.11),plot.margin = margin(1,1,1,1, "cm"),
                                                               legend.title=element_text(size=9),legend.text=element_text(size=7),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1)) 
dev.off()

###############################

library(umap)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/Data.txt",header=T,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("COVID19_NEG", "COVID19_Convalescent", "Hospitalised_Mild","Hospitalised_Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#39540b",COVID19_Convalescent="#6f9035",Hospitalised_Mild="#ccac00",Hospitalised_Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#48690E",COVID19_Convalescent="#8BB443",Hospitalised_Mild="#ffd700",Hospitalised_Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.81, 0.9),plot.margin = margin(1,1,1,1, "cm"),
                                                               legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()


##############################

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/DataFiltered.txt",header=T,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP_Filtered.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP_Filtered.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("COVID19_NEG", "COVID19_Convalescent", "Hospitalised_Mild","Hospitalised_Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/Cluster/UMAP-Filtered.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+stat_ellipse(aes(colour=Group),level=0.95)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#335b76",COVID19_Convalescent="#9eb4cc",Hospitalised_Mild="#ccac00",Hospitalised_Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(COVID19_NEG="#48690E",COVID19_Convalescent="#C6E2FF",Hospitalised_Mild="#ffd700",Hospitalised_Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.81, 0.9),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()

#########
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/MW/Data2Grp.tab",row.names = 1,check.names = FALSE)
T=t(Dat)
write.table(T,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/MW/Data2GrpT.tab",sep="\t",quote = FALSE,col.names = NA)



pval=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/MW/MW_Results_HC_Covid.tab",header = T)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/MW/MW_Results_HC_CovidAdj.tab",sep="\t",col.names = NA, quote = FALSE)

############################ Random Forest

mydata=read.table("/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/Data2GrpT.tab",sep="\t",header=T,check.names = FALSE)

mydata$Group = as.factor(mydata$Group)
library(randomForest)
mtry <- tuneRF(mydata[-1],mydata$Group, ntreeTry=500,stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(best.m)
print(mtry)
set.seed(71)
rf <-randomForest(Group~.,data=mydata, mtry=best.m, importance=TRUE,ntree=200)
print(rf)
imp=importance(rf)
write.table(imp,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/IMP_HC_vs_Covid.txt",sep="\t",col.names = NA, quote = FALSE)


#################

Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/3_grp/Data.txt",row.names = 1,check.names = FALSE)
T=t(Dat)
write.table(T,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/3_grp/DataT.txt",sep="\t",quote = FALSE,col.names = NA)

mydata=read.table("/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/3_grp/DataT.txt",sep="\t",header=T,check.names = FALSE)

mydata$Group = as.factor(mydata$Group)
library(randomForest)
mtry <- tuneRF(mydata[-1],mydata$Group, ntreeTry=500,stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(best.m)
print(mtry)
set.seed(71)
rf <-randomForest(Group~.,data=mydata, mtry=best.m, importance=TRUE,ntree=1000)
print(rf)
imp=importance(rf)
write.table(imp,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/RF/3_grp/IMP.txt",sep="\t",col.names = NA, quote = FALSE)

pred1=predict(rf,type = "prob")
library(ROCR)
perf = prediction(pred1[,2], mydata$Group)


######################## PLS-DA
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Data2Grp.tab",row.names = 1,check.names = FALSE)
Ldat=log2(data)

write.table(Ldat,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Log2data.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Quantile.txt",check.names = FALSE)
M=melt(DAT)
head(M)
Cov=subset(M,Mets="Covid")
hc=subset(M,variable="HC")
?normalizeQuantiles()
ggqqplot(log(Cov$value))

head(data)
Q=normalizeQuantiles(as.matrix(data))
head(Q)

write.table(Q,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Quantile.txt",sep="\t",col.names = NA,quote = FALSE)


library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Info.txt")
group <- paste(sampleinfo$group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.table("/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Limma.txt",sep="\t",row.names=1,header=TRUE,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsCovid=Covid - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCovid",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/sPLS-DA/Result_Limma.txt",sep="\t",quote=FALSE,col.names = NA)


#######################################

library(mixOmics)
data(srbct)
X <- srbct$gene
Y <- srbct$class
summary(Y)

MyResult.splsda <- plsda(X, Y)
plotIndiv(MyResult.splsda)
plotVar(MyResult.splsda)
selectVar(MyResult.splsda, comp=1,)$name 

MyResult.plsda2 <- plsda(X,Y, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 10)
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
list.keepX <- c(5:10,  seq(20, 100, 10))
set.seed(30)
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 3, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)
error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]

############################################# NEW PLAN
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Input.txt",row.names = 1,check.names = FALSE)
L=log2(DAT)
write.table(L,file = "/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Log2Data.txt",sep="\t",col.names = NA,quote = FALSE)

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Metadata.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsMild=Mild - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=Severe - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsConv=Conv - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsConv",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_HCvsConv.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=Severe - Mild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsMild=Mild - Conv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_ConvvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsSevere=Severe - Conv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_ConvvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

##########################################

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Metadata.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsMild=Mild - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=Severe - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=Severe - Mild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


######################## MUVR
library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/X.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("Severe",12))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

######################## HC-Mild

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/XX.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("HC",31))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Mild/Variables.txt",sep="\t",col.names = NA,quote = FALSE)
##################  HC-Severe

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/XX.txt",check.names = FALSE)
Y=c(rep("Severe",12),rep("HC",31))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/HC-Severe/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

######################## Heatmap


############################### Volcano

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = BM),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsMild.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_HCvsMild.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.55,nudge_y = 1.2)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_HCvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_HCvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 1,nudge_y = 0.7)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Limma_MildvsSevere.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Volcano_MildvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = -0.2,nudge_y = 0.24)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/M_S.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/M_S_Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()
#########################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/MetsL.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Temp/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Metadata.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Conv"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))

ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#C6E2FF","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))
colnames(sampleinfo)
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(41, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Heatmap.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

C1=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/C1.txt")
C=melt(C1)
head(C)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Cluster1.pdf")
ggplot(C, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c("HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-3,3)+
  theme(plot.margin = margin(3.5,3.5,3.5,3.5, "cm"),legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_text(size=15,color="black"),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 15,colour = "black"))+
  labs(x="Zscore")
dev.off()

C2=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/C2.txt")
C=melt(C2)
head(C)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/HeatMap/Cluster2.pdf")
ggplot(C, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c("HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-3,3)+
  theme(plot.margin = margin(3.5,3.5,3.5,3.5, "cm"),legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_text(size=15,color="black"),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 15,colour = "black"))+
  labs(x="Zscore")
dev.off()


###################### Correlation

library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Metabolom.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Results.txt",sep="\t",col.names = NA,quote = FALSE)





##################### Bio Marker UMAP

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/BioM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+stat_ellipse(aes(colour=Group),level=0.90,linetype = 2)+
  scale_color_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),axis.text=element_text(size = 12),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()

############# PCA Biomarker
library(PCAtools)
count=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/BioM.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/Metdata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Krushkal/HC_Covd_SigniPCA.txt",row.names = 1,check.names = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/UMAP/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+stat_ellipse(aes(colour=Group),level=0.95)+
  scale_color_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Hospitalized-mild","Hospitalized-severe"),values=c(Mild="#ffd700",Severe="#ffa500"))+
  labs(x="PC1, 61.06%",y="PC2, 12.27%")+theme(axis.title = element_text(size=13),axis.text=element_text(size = 12),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()
?stat_ellipse
##################

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/Sev.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()



data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_IP.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/24_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()



data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_IP.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/Experiment/58_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = "none",plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()


################## Biomarker Proteins MIld-Severe

library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'

XX=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/XX.txt",check.names = FALSE)
Y=c(rep("Mild",29),rep("Severe",12))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

###############


library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Metabolom.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Results.txt",sep="\t",col.names = NA,quote = FALSE)

#################### IPA figures
detach("package:randomForest", unload=TRUE)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_pwy.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_pwy.pdf")
ggplot(data, aes(x=Zscore, y=pval)) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),aes(x=Zscore, y=pval,size=pval),fill="#99aab5",color="#7a8890",pch=21)+
  geom_point(data=subset(data, Zscore <= -2),aes(x=Zscore, y=pval,size=pval),fill="#667BC4",color="#51629c",pch=21)+
  geom_point(data=subset(data, Zscore >=  2),aes(x=Zscore, y=pval,size=pval),fill="#A83131",color="#862727",pch=21)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.3,nudge_x = 0.5,nudge_y = 0.75)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_Toxicity.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Mild_Toxicity.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.3,nudge_x = 0.5,nudge_y = 0.1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Pwy.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Pwy.pdf")
ggplot(data, aes(x=Zscore, y=pval)) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),aes(x=Zscore, y=pval,size=pval),fill="#99aab5",color="#7a8890",pch=21)+
  geom_point(data=subset(data, Zscore <= -2),aes(x=Zscore, y=pval,size=pval),fill="#667BC4",color="#51629c",pch=21)+
  geom_point(data=subset(data, Zscore >=  2),aes(x=Zscore, y=pval,size=pval),fill="#A83131",color="#862727",pch=21)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Toxicity.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/HC_Severe_Toxicity.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 0.2)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/MIld_Severe_disease.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Figure/MIld_Severe_disease.pdf")
ggplot(data, aes(x=Zscore, y=-log10(pval))) + 
  geom_point(data=subset(data, Zscore >-2 & Zscore<2 ),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#99aab5",color="#7a8890")+
  geom_point(data=subset(data, Zscore <= -2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#667BC4",color="#51629c")+
  geom_point(data=subset(data, Zscore >=  2),pch=21,aes(x=Zscore, y=-log10(pval),size=-log10(pval)),fill="#A83131",color="#862727")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,
                   segment.alpha=0.75,box.padding=0.5,nudge_x = 0.5,nudge_y = 0.1)+
  theme(legend.title=element_text(size=12),legend.text=element_text(size=12),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "bottom",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Activation Z-score (IPA)",y="-log10 (P.adj)")+guides(size=guide_legend(override.aes=list(fill="grey",color="grey"),title = "-log10 (P.adj)",nrow = 1))
dev.off()


#########################

library(psych)
library(reshape2)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)

corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Results.txt",sep="\t",col.names = NA,quote = FALSE)
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Melted.txt")
head(DAT)
Re <- acast(DAT,Var1~Var2)
head(Re)
write.table(Re,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.txt",sep="\t",col.names = NA,quote = FALSE)


library(ComplexHeatmap)
library(circlize)
col_corr1= colorRamp2(c(-1, -0.75,-0.5,-0.3, 0.3,0.5,0.75,1), c("#1919ff","#4c4cff","#6f6fff","#9393ff" ,"#ff9999","#ff3232" ,"#ff0000","#cc0000"))
corr=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.txt",check.names = FALSE,row.names = 1)
head(corr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/For_HeatMap.pdf",width = 7,height = 10)
Heatmap(as.matrix(corr),cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width = unit(8, "cm"),show_row_names = TRUE,border = FALSE,
        row_names_gp =gpar(fontsize = 8),height  = unit(10, "cm"),column_names_gp =gpar(fontsize = 8),na_col = "#e6e6e6")
dev.off()

#####################

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Covid_Protein.txt",check.names = FALSE)
mbl=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/MBL.txt",check.names = FALSE)

Res1=corr.test(as.matrix(Met),as.matrix(mbl),use = "pairwise",method="spearman",adjust="none")
Res2=corr.test(as.matrix(Prot),as.matrix(mbl),use = "pairwise",method="spearman",adjust="none")

corr=melt(Res1$r)
pval=melt(Res1$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/MBL_Met.txt",sep="\t",col.names = NA,quote = FALSE)

##########################   Correlation Heatmap
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid_Metabolome.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid_Protein.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")
head(Res$r)
corrplot(Res$r, p.mat = Res$p, sig.level = 0.05,insig = "blank",order="AOE",is.corr = FALSE)
library(corrplot)
warnings()
library(ggcorrplot)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Covid.pdf",width = 10,height = 10)
ggcorrplot(Res$r,p.mat = Res$p,sig.level = 0.1,hc.order = FALSE,insig = "blank",lab_col="black",colors=c("#0000cc","white","#cc0000"),
              tl.cex=5,tl.srt=90,legend.title="Correlation",outline.color="#e6e6e6")
dev.off()
head(DD1)
write.table(DD1$data,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.txt",sep="\t",col.names = NA,quote = FALSE)

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.txt")
head(DAT)
library(reshape2)
Re <- acast(DAT,factor(Var1,levels=unique(Var1))~factor(Var2,levels=unique(Var2)))

head(Re)
library(ComplexHeatmap)
library(circlize)
col_corr= colorRamp2(c(-1, -0.75,-0.5,-0.3,0, 0.3,0.5,0.75,1), c("#1919ff","#4c4cff","#6f6fff","#9393ff" ,"white","#ff9999","#ff3232" ,"#ff0000","#cc0000"))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/Correlation.pdf",width = 7,height = 10)
Heatmap(as.matrix(Re),col=col_corr,cluster_rows=FALSE,cluster_columns = FALSE,name="Spearman",width  = unit(8, "cm"),show_row_names = TRUE,border = FALSE,
        row_names_gp =gpar(fontsize = 8),height  = unit(23, "cm"),column_names_gp =gpar(fontsize = 8),na_col = "#e6e6e6")
dev.off()

############ For soham

Ol=read.delim("/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/RNA.txt",check.names=FALSE)
Ref=read.delim("/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/parp.txt",check.names = FALSE)
Res=corr.test(as.matrix(Ol),as.matrix(Ref),use = "pairwise",method="spearman",adjust="none")
head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/FinalFigures/Gender_Deseq2/For_Soham/res_RNA.txt",sep="\t",col.names = NA,quote = FALSE)

################  Sankey Plot Mild-Severe

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/pwy.txt")
library(ggplot2)
library(ggalluvial)
head(XX)
geom_stratum(fill = c("#975d72","#2675a0","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2"),width = 0.2) +
  
col=c(rep("black",24),rep("#92bacf",28))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/IPA/Pathway.pdf",height = 10,width = 7.5)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Metabolites)) +
  geom_alluvium(aes(),fill="#ec8e56")+#scale_fill_manual(values = c("#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56"))+
  geom_stratum(fill = c("#98b3c2","#98b3c2","#98b3c2","#98b3c2","#98b3c2",
                        "#98b3c2","#98b3c2","#98b3c2","#98b3c2","#98b3c2",
                        "#98b3c2","#98b3c2","#98b3c2","#98b3c2","#A83131",
                        "#667BC4","#667BC4","#667BC4","#A83131","#A83131",
                        "#A83131","#667BC4","#A83131"),width = 0.2) +
  geom_label(stat = "stratum", infer.label = TRUE,size=4) +
  scale_x_discrete(limits = c("Pathway", "Metabolites"), expand = c(.3, .3)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank())+theme(legend.position = "none",plot.margin = margin(3,0.5,3,0.5, "cm"),
                                                                          axis.text.x = element_text(size=15,color="black"))
dev.off()

############# Violin
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/BioMarker.txt",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Mets,y=`eicosanedioate (C20-DC)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "eicosanedioate (C20-DC)")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P2=ggplot(Bean,aes(x=Mets,y=mannose,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "mannose")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2.5, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P3=ggplot(Bean,aes(x=Mets,y=`eicosenedioate (C20:1-DC)*`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "eicosenedioate (C20:1-DC)*")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 3, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P4=ggplot(Bean,aes(x=Mets,y=`hydantoin-5-propionate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "hydantoin-5-propionate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 3, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P5=ggplot(Bean,aes(x=Mets,y=`6-bromotryptophan`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-bromotryptophan")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 2, by = 1.2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

#P6=ggplot(Bean,aes(x=Mets,y=`1-palmitoyl-2-oleoyl-GPC (16:0||18:1)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "1-palmitoyl-2-oleoyl-GPC")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1.5, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P7=ggplot(Bean,aes(x=Mets,y=`4-hydroxyphenylacetate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "4-hydroxyphenylacetate")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 4, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

P8=ggplot(Bean,aes(x=Mets,y=`6-oxopiperidine-2-carboxylate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-oxopiperidine-2-carboxylate")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 3, by = 1.2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position =c(1.7,0.5),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 13),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))

#P6=ggplot(Bean,aes(x=Mets,y=`1-palmitoyl-2-oleoyl-GPC (16:0||18:1)`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "1-palmitoyl-2-oleoyl-GPC")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1.5, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(1.4,1.1),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 13))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/BioMarker.pdf",height = 10,width =8)
ggarrange(P2,P1,P5,P3,P7,P4,P8,nrow = 4,ncol=2)+theme(plot.margin = margin(0.5,7,0.5,0.5, "cm"))
dev.off()

############### Violin pathway

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/pathway.txt",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Mets,y=`alpha-ketobutyrate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "alpha-ketobutyrate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 2.5, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=20),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))




P3=ggplot(Bean,aes(x=Mets,y=`glucose`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "D-glucose")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 2, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P4=ggplot(Bean,aes(x=Mets,y=`formiminoglutamate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "formiminoglutamic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-2.5, 2.5, by = 1.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))


P5=ggplot(Bean,aes(x=Mets,y=`hydantoin-5-propionate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "hydantoin-5-propionate")+
  theme_bw()+scale_y_continuous(breaks = seq(-3.5, 3, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P6=ggplot(Bean,aes(x=Mets,y=`indoleacetate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "indoleacetic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 5, by = 2))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P7=ggplot(Bean,aes(x=Mets,y=`alanine`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "L-alanine")+
  theme_bw()+scale_y_continuous(breaks = seq(-1.5, 1, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P8=ggplot(Bean,aes(x=Mets,y=`tryptophan`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "L-tryptophan")+
  theme_bw()+scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.6))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P9=ggplot(Bean,aes(x=Mets,y=`quinolinic`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "quinolinic acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-2, 2.5, by = 1))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))

P2=ggplot(Bean,aes(x=Mets,y=`citrate`,fill=Mets,color=Mets))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "citric acid")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 1, by = 0.5))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#324909",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position ="none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=18),plot.title = element_text(hjust = 0.5,size = 20))
library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Pathway.pdf",height = 10,width =30)
ggarrange(P1,P4,P5,P3,P9,P2,P6,P7,P8,nrow = 1)+theme(plot.margin = margin(5,0,5,0, "cm"))
dev.off()


#################################


Bean=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Data.txt",check.names = FALSE)


P1=ggplot(Bean,aes(x=Group,y=HGF,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+
  labs(y="NPX",title = "HGF")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=9),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P2=ggplot(Bean,aes(x=Group,y=PTN,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "PTN")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 8, by = 4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P3=ggplot(Bean,aes(x=Group,y=CXCL13,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL13")+
  theme_bw()+scale_y_continuous(breaks = seq(6, 14, by = 2))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P4=ggplot(Bean,aes(x=Group,y=`MCP-3`,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MCP-3")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 10, by = 3))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")

P5=ggplot(Bean,aes(x=Group,y=IL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "IL12")+
  theme_bw()+scale_y_continuous(breaks = seq(4, 12, by = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),legend.position = "none")




P9=ggplot(Bean,aes(x=Group,y=CXCL12,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "CXCL12")+
  theme_bw()+scale_y_continuous(breaks = seq(0, 2, by = 0.4))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(color="black",size=9),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.text.x = element_text(size = 9,color="black"),
        legend.position = c(1.6,1.3))




library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Violin.pdf",height = 5,width =9)
ggarrange(P1,P2,P3,P4,P5,P9,nrow = 2,ncol = 3)+theme(plot.margin = margin(2,5,2,2, "cm"))
dev.off()

##############################################

ip=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Rank.txt")
head(ip)
ip$Mets <- factor(ip$Mets, levels = ip$Mets)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Rank.pdf")
ggplot(ip, aes(y=Mets)) + 
  geom_point(data=ip,aes(x=1,y=Mets,size=rank,color=rank))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#191978",high="#7f7fb4",breaks=c(1,10,50,200,400))+
  scale_y_discrete(position = "right")+theme_bw()+coord_fixed(ratio = 0.95)+scale_size(range = c(5,8),breaks=c(1,10,50,200,400))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(5,2,5.4,2.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(3, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Rank",nrow = 1),color=guide_legend(title = "Rank",nrow=1))
dev.off()



ip=read.delim("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Rank.txt")
head(ip)
ip$Mets <- factor(ip$Mets, levels = ip$Mets)
pdf("/home/anoop/Desktop/COVID_Omics/New_Analysis/MUVR/Rank.pdf")
ggplot(ip, aes(y=Mets)) + 
  geom_point(data=ip,aes(x=1,y=Mets,size=rank,color=rank))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#65880f",high="#a2b76f",breaks=c(1,5,10,20,50))+
  scale_y_discrete(position = "right")+theme_bw()+coord_fixed(ratio = 0.99)+scale_size(range = c(5,8),breaks=c(1,5,10,20,50))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(5,3.5,7,2.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(1.8, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Rank",nrow = 1),color=guide_legend(title = "Rank",nrow=1))
dev.off()

#################### Donut

data <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Donut.txt")
head(data)
data$Level <- factor(data$Level, levels = data$Level)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Donut.pdf",width = 10,height = 10)
ggdonutchart(data, size=0.7,"Percent", label = "Label",fill="Category",color = "white")+theme(plot.margin = margin(4,10,4,0, "cm"),axis.text = element_text(size=22),
                                                                                              legend.position = c(1.35,0.5),legend.text = element_text(size = 15),
                                                                                              legend.title = element_blank())
dev.off()


#####
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/Input.txt")
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MBL",y="ng/mL")+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = c(1.8, 0.5))+scale_y_continuous(breaks = seq(0.5, 12, by =3))

P2=ggplot(Bean,aes(x=Group,y=Neopterin,fill=Group,color=Group))+geom_violin()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "Neopterin",y="nmol/L")+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = c(1.4, 0.5))+scale_y_continuous(breaks = seq(4, 55, by =15))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Maike/Violin.pdf",height = 5,width =10)
ggarrange(P1,P2,nrow = 1)+theme(plot.margin = margin(2,6,2,0.5, "cm"))
dev.off()

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/MBL.txt")
pdf("/home/anoop/Desktop/COVID_Omics/Maike/MBL_newActualValue.pdf",height = 5,width =10)
ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",mild="#ccac00",severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",mild="#ffd700",severe="#ffa500"))+labs(title = "MBL",y = expression(mu*"g/mL"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),legend.text = element_text(size = 9),
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),plot.margin = margin(2,13,2,0.5, "cm"),
        panel.border = element_blank(),legend.position = c(1.25, 0.5))+scale_y_continuous(breaks = seq(0, 5, by =1))
dev.off()


####################### NEW Limma with Co-factors

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/MetaData.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsConv=GroupConv - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsConv",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsConv.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsMild=GroupMild - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsSevere=GroupSevere - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(ConvvsMild=GroupMild - GroupConv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_CONVvsMild.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(ConvvsSevere=GroupSevere - GroupConv,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ConvvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_CONVvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(MildvsSevere=GroupSevere - GroupMild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

################################
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/MetaDataHC_Covid.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsCovid=GroupCovid - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCovid",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Limma_HCvsCovid.txt",sep="\t",quote=FALSE,col.names = NA)

########################  3 Group with co-variates

sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/MetaData.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsMild=GroupMild - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsMild",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsMild.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(HCvsSevere=GroupSevere - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(MildvsSevere=GroupSevere - GroupMild,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="MildvsSevere",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_MildvsSevere.txt",sep="\t",quote=FALSE,col.names = NA)


############# Volcano plots NEW

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsMild2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_HCvsMild.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = -0.2,nudge_y = 0.24)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_HCvsSevere2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_HCvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = 0.1,nudge_y = 0.15)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Limma_MildvsSevere2.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/3_Group/Volcano_MildvsSevere.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_label_repel(aes(label = BM),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0.2,nudge_x = 0.2,nudge_y = 0.1)+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/MetaDataHC_Covid.txt")

design <- model.matrix(~0+Age+BMI+Gender+Group,data = sampleinfo)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
head(design)
cont.matrix <- makeContrasts(HCvsCovid=GroupCovid - GroupHC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCovid",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Adjusted/Limma_HCvsCovid.txt",sep="\t",quote=FALSE,col.names = NA)


####################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/MetsL.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Heatmap1.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Design.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Conv"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))

head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#8BB443","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))

MET=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Mets.txt",row.names = 1)
ha2 = rowAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Super_Pathway = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Super_Pathway=c("Amino Acid"="#ed8d18","Carbohydrate"="#017399","Cofactors and Vitamins"="#f2af5d","Energy"="#02ace5",
                                                  "Lipid"="#8e540e","Nucleotide"="#e5dfa4","Peptide"="#7f7c5b")))


ha3 = rowAnnotation(foo = anno_mark(at = c(38,55,100,119,123,150,32,56,58,61,47,111,37,101,122,152,50),
                                    labels_gp = gpar(fontsize=20),lines_gp = gpar(col="black"),link_height = unit(35, "mm"),link_width=unit(12, "mm"),
                                    labels = c("aspartate","glutamate","phenylalanine","fructose","mannose","alpha-ketoglutarate","alanine","glutamine",
                                               "glycine","histidine","cysteine","tryptophan","arginine","proline","glycerate","citrate","cystine")))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           right_annotation = ha3,left_annotation= ha2,top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 4),height  = unit(41, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",27),rep("d_Covid_severe",12)))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/Heatmap.pdf",height = 20,width =25)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

packageVersion("GGally")
Cls=row_order(H1)
row.names(Zscore)

??annotation_legend_param
clu = t(t(row.names(Zscore[row_order(H1)[[2]],])))
write.table(as.data.frame(clu),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Heatmap/CL2.txt",sep="\t",quote = FALSE,col.names = NA)

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/superCnt.txt")

library(ggplot2)
data$Super <- factor(data$Super, levels = data$Super)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/SuperPathway.pdf")
ggplot(data, aes(x = Type,y = Perc, fill = Super)) +theme_bw()+labs(y="Percetage")+guides(fill=guide_legend(title="Super pathway"))+
  geom_bar(stat="identity")+scale_fill_manual(values = c("Lipid"="#8e540e","Amino Acid"="#ed8d18","Cofactors and Vitamins"="#f2af5d","Peptide"="#7f7c5b",
                                                         "Nucleotide"="#e5dfa4","Carbohydrate"="#017399","Energy"="#02ace5","Partially Characterized Molecules"="#35ccff"))+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(2.5,3,2.5,3, "cm"))


dev.off()

lipid=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/Lipid.txt")
head(lipid)
lipid$Type <- factor(lipid$Type, levels = lipid$Type)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Bargraph/Lipid.pdf")
ggplot(lipid, aes(x = Group,y = Perc, fill = Type)) +theme_bw()+labs(y="Percetage")+guides(fill=guide_legend(title="Sub pathway"))+
  geom_bar(stat="identity")+scale_fill_manual(values = c("Fatty Acid, Dicarboxylate"="#994500","Lysophospholipid"="#e56800","Androgenic Steroids"="#ff8f32",
                                                         "Sphingomyelins"="#ffb97f","Monoacylglycerol"="#665600","Diacylglycerol"="#998100","Plasmalogen"="#ccac00",
                                                         "Phosphatidylcholine (PC)"="#fdd017",
                                                         "Long Chain Polyunsaturated Fatty Acid (n3 and n6)"="#746776","Phosphatidylethanolamine (PE)"="#baa5bc",
                                                         "Hexosylceramides (HCER)"="#e9cfec","Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)"="#4e599e",
                                                         "Secondary Bile Acid Metabolism"="#6473cb","Fatty Acid, Monohydroxy"="#8c99e7",
                                                         "Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)"="#0aa091",
                                                         "Fatty Acid Metabolism (Acyl Choline)"="#25cebd",
                                                         "Other"="#9ee9e1"))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),plot.margin = margin(2.5,2,2.5,1.5, "cm"))
dev.off()


############# Sankey  HC-Covid After adjusted

XX=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid.txt")
library(ggplot2)
library(ggalluvial)
head(XX)

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid_pwy.pdf",height = 10,width = 10)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Metabolite)) +
  geom_alluvium(aes(),fill="#ec8e56")+#scale_fill_manual(values = c("#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56","#ec8e56"))+
  geom_stratum(width = 0.1) +
  geom_label(stat = "stratum", infer.label = TRUE,size=2) +
  scale_x_discrete(limits = c("Pathway", "Metabolites"), expand = c(.5, .5)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),panel.grid = element_blank())+theme(legend.position = "none",plot.margin = margin(1,0.5,1,0.5, "cm"),
                                                                                                       axis.text.x = element_text(size=15,color="black"))
dev.off()



################## NEW Experiment Soham

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Caco2_Meta.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)


ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("UnTreated"="#48690E","Treated"="orange")))
colnames(sampleinfo)
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(10, "cm"),width  = unit(10, "cm"),
           column_split =c(rep("a_UT",3),rep("b_T",3)))
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/Viral_Caco2.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

############ ISG Volcano plot
library(ggrepel)
data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Caco2_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Caco2_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-1.5, 2.5)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Calu3_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/ISG/Calu3_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-2, 2)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/293T_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/293T_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-1, 1)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/Huh7_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/293T_Huh7/ISG/Huh7_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + xlim(-0.5, 0.5)+
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "#c2f542", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()



ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/NonMetabolic.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/NonMetabolic.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#191978",high="#7f7fb4",breaks=c(1E-20,1E-10,0.000001,0.0001,0.001,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2,6,4,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/Metabolic.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Calu3/Metabolic.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#664c7f",high="#ead6ff",breaks=c(1E-20,1E-10,0.000001,0.0001,0.001,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(0,6.2,0,5.7, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(15,0.5),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Caco2/Input.txt")
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/NewExperiment/Calu3_Caco2/KEGG_Caco2/Caco2.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Ratio,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#05a167",high="#9bd9c2",breaks=c(0.0001,0.001,0.05,0.1,0.3,0.5))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,6,5,4.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(12.5,0.5),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Ratio(%)",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()

########## pairs plot

library(ggplot2)
library(GGally)
data <- data.frame( var1 = 1:100 + rnorm(100,sd=20), v2 = 1:100 + rnorm(100,sd=27), v3 = rep(1, 100) + rnorm(100, sd = 1)) 
data$v4 = data$var1 ** 2 
data$v5 = -(data$var1 ** 2)
?ggpairs
both=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/ggpair/Data.txt")
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/Both.pdf")
ggpairs(both, title="correlogram with ggpairs()") 
dev.off()
head(flea)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/ggpair/Corr.pdf")
ggpairs(both, mapping = aes(fill=Group,color = Group,alpha = 0.5),showStrips=NULL,lower = list(continuous = "smooth",size=5),
        upper = list(continuous = wrap("cor", size = 3)),diag = list(continuous = "barDiag")
        )+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(size=9))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))
dev.off()

####### Correlation Biomarkers 3 groups

library(psych)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_HC.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_HC.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_HC_1.txt",sep="\t",col.names = NA,quote = FALSE)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_Mild.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_Mild.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Mild_1.txt",sep="\t",col.names = NA,quote = FALSE)

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Meta_Severe.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Prot_Severe.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Severe_1.txt",sep="\t",col.names = NA,quote = FALSE)


Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Severe.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Severe_MBL.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Severe_MBL.txt",sep="\t",col.names = NA,quote = FALSE)

###

Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Covid.txt",check.names=FALSE)
Prot=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Covid_MBL.txt",check.names = FALSE)
Res=corr.test(as.matrix(Met),as.matrix(Prot),use = "pairwise",method="spearman",adjust="none")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value
?corr.test
write.table(corr,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/BIomarker/FINAL/Results_Covid_MBL_1.txt",sep="\t",col.names = NA,quote = FALSE)


################## New violin
library(ggplot2)
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/Maike/NewViolin.txt")
head(Bean)

P1=ggplot(Bean,aes(x=Group,y=Glucose,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glucose")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-4, 2, by =0.65))

P2=ggplot(Bean,aes(x=Group,y=Mannose,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Mannose")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1))


P3=ggplot(Bean,aes(x=Group,y=MBL,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "MBL",y = expression(mu*"g/mL"))+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(0, 5, by =1.25))

P4=ggplot(Bean,aes(x=Group,y=Alanine,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Alanine")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.65))


P5=ggplot(Bean,aes(x=Group,y=Tryptophan,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Tryptophan")+
  theme_bw()+stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.6))


P6=ggplot(Bean,aes(x=Group,y=Glutamate,fill=Group,color=Group))+geom_boxplot()+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamate")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =1))

my_comparisons = list( c("Mild", "Severe"))
?stat_compare_means
library(ggsignif)
library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Maike/NewViolin.pdf",height = 5,width =9)
ggarrange(P1,P2,P3,P4,P5,P6,nrow = 2,ncol = 3)+theme(plot.margin = margin(1.5,5,1.5,2, "cm"))
dev.off()


#############

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/Input.txt",header = TRUE)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Limma_3_grp/NEW/Sankey/HC_Covid.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=pval,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#fce79a",high="#b3962d")+
  scale_y_discrete(position = "right")+theme_bw()+scale_size(range = c(5,7))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position = c(6, -0.1),
        legend.box="vertical",legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),
        panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(4.5,5,4.5,5.2, "cm"))+
  guides(color=guide_legend(nrow =1 ,title = "-log10(Adj.Pvalue)"), size = guide_legend(nrow = 2,title = "-log10(Adj.Pvalue)"))
dev.off()
################## Supplement Box
Bean=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Suppl.txt",header = TRUE)  
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=glycine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glycine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 9),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 2, by =0.75))


P2=ggplot(Bean,aes(x=Group,y=proline,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Proline",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.65))

P3=ggplot(Bean,aes(x=Group,y=tryptophan,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Tryptophan",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.55))


P4=ggplot(Bean,aes(x=Group,y=alanine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Alanine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =0.55))


P5=ggplot(Bean,aes(x=Group,y=histidine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Histidine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.4))


P6=ggplot(Bean,aes(x=Group,y=glutamine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.35))


P7=ggplot(Bean,aes(x=Group,y=arginine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Arginine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.75))





P9=ggplot(Bean,aes(x=Group,y=aspartate,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Aspartate",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 2, by =0.75))

P10=ggplot(Bean,aes(x=Group,y=cystine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Cystine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.5))

P11=ggplot(Bean,aes(x=Group,y=phenylalanine,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Phenylalanine",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 1, by =0.45))


P12=ggplot(Bean,aes(x=Group,y=fructose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Fructose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1.25))


P13=ggplot(Bean,aes(x=Group,y=mannose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Mannose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 3, by =1))


P14=ggplot(Bean,aes(x=Group,y=glucose,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glucose",y="Log2(Measurement)")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-1, 2, by =0.65))

P8=ggplot(Bean,aes(x=Group,y=glutamate,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",HC_cov="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",HC_cov="#8BB443",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Glutamate",y="Log2(Measurement)")+
  theme_bw()+guides(color=guide_legend(nrow =1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 9),legend.position = c(2.5,-0.2),legend.title=element_blank(),
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())+scale_y_continuous(breaks = seq(-2, 1, by =0.85))
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Violin/Suppl.pdf",height = 5,width =12)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,nrow = 2,ncol = 7)+theme(plot.margin = margin(1,0,1.5,0, "cm"))
dev.off()

############## Facs
library(reshape2)
C1=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",header = TRUE)
head(C1)
C=melt(C1)
library(ggplot2)
library(ggridges)

D1=subset(C,variable=="CD8_GLUT")
P1=ggplot(D1, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(90,102), breaks = seq(90,102, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "CD8 GLUT")
?stat_density_ridges
D2=subset(C,variable=="IM_GLUT")
P2=ggplot(D2, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(-2,9), breaks = seq(-2,9, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "IM GLUT")

D3=subset(C,variable=="NCM_GLUT")
P3=ggplot(D3, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(-2,9), breaks = seq(-2,9, by = 5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "NCM GLUT")

D4=subset(C,variable=="NCM_XCT")
P4=ggplot(D4, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(0,120), breaks = seq(0,120, by = 50))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "NCM XCT")

D5=subset(C,variable=="CM_XCT")
P5=ggplot(D5, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(99,100.2), breaks = seq(99,100.2, by = 0.5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "CM XCT")

D6=subset(C,variable=="IM_XCT")
P6=ggplot(D6, aes(x=value,y=Mets,fill=Mets,alpha=0.1,color=Mets))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Mets),quantile_lines = TRUE,quantiles = 0.5)+scale_x_continuous(limits = c(99,100.2), breaks = seq(99,100.2, by = 0.5))+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"),plot.margin = margin(0.1,2,0.1,2, "cm"))+
  labs(x="Zscore",title = "IM XCT")

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/ubha/Facs.pdf",width = 10,height = 10)
ggarrange(P1,P4,P2,P5,P3,P6,nrow = 3,ncol = 2)+theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()



C1=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",header = TRUE)
data=read.table("/home/anoop/Desktop/COVID_Omics/ubha/Data.txt",sep="\t",header=TRUE)
rnames <- data[,1] 
mat_data <- data.matrix(data[,2:ncol(data)])
S=scale(mat_data, center = TRUE, scale = TRUE)
rownames(S) <- rnames

C=melt(S)
library(ggplot2)
library(ggridges)
head(C)
D1=subset(C,Var2=="CD8_GLUT")
P1=ggplot(D1, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "bottom",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "CD8 GLUT")

D2=subset(C,Var2=="IM_GLUT")
P2=ggplot(D2, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "IM GLUT")

D3=subset(C,Var2=="NCM_GLUT")
P3=ggplot(D3, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "NCM GLUT")


D4=subset(C,Var2=="NCM_XCT")
P4=ggplot(D4, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "NCM XCT")

D5=subset(C,Var2=="CM_XCT")
P5=ggplot(D5, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "CM XCT")


D6=subset(C,Var2=="IM_XCT")
P6=ggplot(D6, aes(x=value,y=Var1,fill=Var1,alpha=0.1,color=Var1))+scale_fill_manual(values = c(FMT="grey",HC_cov="#8BB443","HC"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c(FMT="grey",HC_cov="#6f9035","HC"="#9eb4cc","Mild"="#ccac00","Severe"="#cc8400"))+
  stat_density_ridges(aes(fill = Var1),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+
  theme(legend.position = "none",axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(size = 10,colour = "black"))+
  labs(x="Zscore",title = "IM XCT")


pdf("/home/anoop/Desktop/COVID_Omics/ubha/FacsScaled.pdf",width = )
ggarrange(P1,P4,P2,P5,P3,P6,nrow = 3,ncol = 2)+theme(plot.margin = margin(1,1,1,1, "cm"))
dev.off()


###########################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Input.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Design.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)


head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("Healthy Control (HC)"="#48690E","HC (CoV-2 Ab+)"="#8BB443","Hospitalized-mild"="#ffd700","Hospitalized-severe"="#ffa500")))

MET=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Pwy.txt",row.names = 1)
ha2 = rowAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Pwy = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                    col = list(Pwy=c("Glycolysis, Gluconeogenesis, and Pyruvate Metabolism"="#0284d0","Fructose, Mannose and Galactose Metabolism"="#f07aaa","TCA Cycle"="#045463")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           left_annotation= ha2,top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(18, "cm"),width  = unit(40, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",29),rep("d_Covid_severe",12)))


col_fun_lfc = colorRamp2(c(-1, -0.05,0, 0.5,1), c("#065535","#508871" ,"white", "#ff7f7f","#e50000"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/LFC.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(2, "cm"),show_row_names = TRUE,
           heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           row_names_gp =gpar(fontsize = 20),height  = unit(18, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 
t=H1+H2


pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/Heatmap.pdf",height = 20,width =25)
draw(t,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

####### Dotplot

data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.txt",header = TRUE,check.names = FALSE)
head(data)
P1=ggplot(data, aes(x=Group, y=Fructose,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Fructose")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits=c(-2,4),breaks = seq(-2, 4, by =1.75))


P2=ggplot(data, aes(x=Group, y=Lactate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Lactate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =1))


P3=ggplot(data, aes(x=Group, y=Pyruvate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Pyruvate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-0.9, 1.5),breaks = seq(-1, 1, by =.65))


P4=ggplot(data, aes(x=Group, y=Citrate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Citrate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.5),breaks = seq(-1, 1, by =.65))

P5=ggplot(data, aes(x=Group, y=`Aconitate [cis or trans]`,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Aconitate [cis or trans]")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =.8))


P6=ggplot(data, aes(x=Group, y=`Alpha-ketoglutarate`,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = expression(alpha*"-ketoglutarate"))+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 2.5),breaks = seq(-2, 3, by =1.3))

P7=ggplot(data, aes(x=Group, y=Fumarate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Fumarate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 1.6),breaks = seq(-5, 1, by =1.7))


P8=ggplot(data, aes(x=Group, y=Malate,fill=Group,color=Group)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Malate")+
  stat_summary(fun = "median", colour = "blue", size = 2, geom = "point")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.2),breaks = seq(-1, 1, by =.55))


P9=ggplot(data, aes(x=Group, y=Succinate,fill=Group,color=Group,size=6)) + theme_bw()+geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="log2 Measurement",title = "Succinate")+
  theme(legend.position = c(1.5,0.5),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 2),breaks = seq(-5, 1.5, by =2))+
  stat_summary(fun = "median", color="blue", fill = "blue", size = 2, geom = "point")

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.pdf",width = 13,height = 10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,nrow = 2,ncol = 5)+theme(plot.margin = margin(6,0,6,0, "cm"))
dev.off()

###########


data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/DotPlot.txt",header = TRUE,check.names = FALSE)
head(data)
P1=ggplot(data, aes(x=Group, y=Fructose,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Fructose")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits=c(-2,4),breaks = seq(-2, 4, by =1.75))


P2=ggplot(data, aes(x=Group, y=Lactate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Lactate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =1))


P3=ggplot(data, aes(x=Group, y=Pyruvate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Pyruvate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-0.9, 1.5),breaks = seq(-1, 1, by =.65))


P4=ggplot(data, aes(x=Group, y=Citrate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Citrate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.5),breaks = seq(-1, 1, by =.65))

P5=ggplot(data, aes(x=Group, y=`Aconitate [cis or trans]`,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Aconitate [cis or trans]")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 1.5),breaks = seq(-2, 1, by =.8))


P6=ggplot(data, aes(x=Group, y=`Alpha-ketoglutarate`,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = expression(alpha*"-ketoglutarate"))+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-2, 2.5),breaks = seq(-2, 3, by =1.3))

P7=ggplot(data, aes(x=Group, y=Fumarate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Fumarate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 1.6),breaks = seq(-5, 1, by =1.7))


P8=ggplot(data, aes(x=Group, y=Malate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Malate")+
  theme(legend.position = "none",axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-1, 1.2),breaks = seq(-1, 1, by =.55))


P9=ggplot(data, aes(x=Group, y=Succinate,fill=Group,color=Group)) + theme_bw()+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  labs(y="log2 Measurement",title = "Succinate")+
  theme(legend.position = c(1.5,0.5),axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 10),axis.title.y = element_blank(),
        panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size = 10,colour = "black"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#8BB443","HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500"))+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values = c("HC_conv"="#6f9035","HC"="#39540b","Mild"="#ccac00","Severe"="#cc8400"))+
  scale_y_continuous(limits = c(-5, 2),breaks = seq(-5, 1.5, by =2))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Heatmap_Glycolysis/BoxPlot.pdf",width = 13,height = 10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,nrow = 2,ncol = 5)+theme(plot.margin = margin(6,0,6,0, "cm"))
dev.off()


Data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/Significant_Result.txt",header = TRUE,row.names = 1)
head(Data)
Pval=Data$padj
myPval <- Data["padj"]
myLFC<-Data["log2FoldChange"]

head(myLFC)
BiocManager::install("rsbml")
n
library(libSBML)
library(piano)
install.packages("/home/anoop/Tools/libSBML_5.19.0.tar.gz", repos = NULL, type="source")
library(rsbml)
myGsc <- loadGSC("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Human-GEM.xml")


library(psych)
library(reshape2)
Met=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Correlation/ForShuba/Severe.txt",check.names=FALSE)

Res=corr.test(as.matrix(Met),use = "pairwise",method="spearman",adjust="none")
melt(Res$r)
melt(Res$p)

##################### ISG heatmap for soham

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG_LCPM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG_LCPM_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG_LCPM_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Figures/Metadata.txt",row.names = 1,header = TRUE)

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Convalescent"="#C6E2FF","hosp_mild"="#ffd700","hosp_severe"="#ffa500"))

head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("HC"="#48690E","Convalescent"="#8BB443","hosp_mild"="#ffd700","hosp_severe"="#ffa500")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,row_split = (MET$Category),cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           ,top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(45, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("a_HC",21),rep("b_Antibody_positive",10),rep("c_Covid_mild",26),rep("d_Covid_severe",11)))

cmp=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/comp.txt",row.names = 1,check.names = FALSE)
ha3 = HeatmapAnnotation(df = cmp,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Comparison = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                        col = list(Comparison=c("HC-vs-hosp_mild"="#b2b2b2","HC-vs-hosp_severe"="#808080","hosp_mild-vs-hosp_severe"="#595959")))

col_fun_lfc = colorRamp2(c(2, 1, 0, -1, -2), c("#98120c","#fe1f14" ,"#FFFCC9","#6666ff","#4747b2"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG_LFC.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(2, "cm"),show_row_names = TRUE,show_column_names = FALSE,
           top_annotation = ha3,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           row_names_gp =gpar(fontsize = 5),height  = unit(45, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 
t=H1+H2


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/ISG/ISG.pdf",height = 20,width =17)
draw(t,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


#####################################
library(package)
library(ggplot2)
data=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Mild_TPM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Mild_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Mild_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Mild_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.2,0.9),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Severe_TPM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Severe_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Severe_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Severe_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.81,0.9),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()


library(ggrepel)
data=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/HC_TPM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/HC_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/HC_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/HC_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.2,0.1),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()


library(ggrepel)
data=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Conv_TPM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Conv_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Conv_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Conv_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#9eb4cc",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#C6E2FF",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.2,0.1),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()


#######################################

ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C1/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C1/C1.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#39540b",high="#91a56e")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=7,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(2,6.2,4,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C2/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C2/C2.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#7b6121",high="#b8a67c",breaks=c(1E-8,0.003,0.004,0.001,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,5,8,6.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.3),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C3/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C3/C3.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#cc8400",high="#ffc04c")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,6.2,3,4, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4, -0.12),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C4/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C4/C4.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#e5c100",high="#ffe34c")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,5,8,5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(3, -0.3),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C5/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C5/C5.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#cc8400",high="#ffc04c")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,5,3,4, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(3, -0.12),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()


ip=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C6/Input.txt",header = TRUE)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/KEGG_C6/C6.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Overlap,color=Padj))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#cc8400",high="#ffc04c")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(1.8,4,3,4.8, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4, -0.12),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()

####################

pval=read.delim("/home/anoop/Desktop/MonoMDC_Panel/MW_Result.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P_value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/MonoMDC_Panel/MW_Result_Adjusted.txt",sep="\t",col.names = NA, quote = FALSE)


###################################### Covid Network Fusion New patient group DEseq2 ##########

library(DESeq2)
count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus2/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls1",22),rep("Cls2",19)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
?results
res=results(ds,c("condition","Cls2","Cls1"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus2/Results.txt",sep="\t", quote=FALSE,col.names = NA)

count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus3/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls1",22),rep("Cls3",10)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Cls3","Cls1"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus3/Results.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus4/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls1",22),rep("Cls4",17)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Cls4","Cls1"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus1_vs_Clus4/Results.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus2_vs_Clus3/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls2",19),rep("Cls3",10)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Cls3","Cls2"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus2_vs_Clus3/Results.txt",sep="\t", quote=FALSE,col.names = NA)

count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus2_vs_Clus4/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls2",19),rep("Cls4",17)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Cls4","Cls2"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus2_vs_Clus4/Results.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus3_vs_Clus4/Input.txt",row.names = 1, check.names = FALSE)
samples=data.frame(condition=c(rep("Cls3",10),rep("Cls4",17)))
ds=DESeqDataSetFromMatrix(countData = count,colData = samples,design = ~condition)
ds=DESeq(ds)
res=results(ds,c("condition","Cls4","Cls3"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/DGE/Clus3_vs_Clus4/Results.txt",sep="\t", quote=FALSE,col.names = NA)


########################################################### Metabolomics

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/MetaData.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(Cls1vsCls2=Cls2 - Cls1,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls1vsCls2",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls1vsCls2.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(Cls1vsCls3=Cls3 - Cls1,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls1vsCls3",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls1vsCls3.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(Cls1vsCls4=Cls4 - Cls1,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls1vsCls4",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls1vsCls4.txt",sep="\t",quote=FALSE,col.names = NA)


cont.matrix <- makeContrasts(Cls2vsCls3=Cls3 - Cls2,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls2vsCls3",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls2vsCls3.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(Cls2vsCls4=Cls4 - Cls2,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls2vsCls4",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls2vsCls4.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(Cls3vsCls4=Cls4 - Cls3,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="Cls3vsCls4",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Metabolome/Cls3vsCls4.txt",sep="\t",quote=FALSE,col.names = NA)

###########################

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Healthy_Covid/Metadata.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Healthy_Covid/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(HCvsCov=Covid - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsCov",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/Healthy_Covid/Results.txt",sep="\t",quote=FALSE,col.names = NA)




################################### NEW FIGURES

#####################################

library(ggplot2)
library(umap)

data=read.delim("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Fitered_TPM.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/TPM_Filtered_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/TPM_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Cohort <- factor(dat$Cohort, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/RNASeq_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.8,0.1),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()


dat=read.delim("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Coding_LCPM_UMAP.txt",row.names = 1,check.names = FALSE)
dat$Cohort <- factor(dat$Cohort, levels=c("HC", "Conv", "Mild","Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Coding_LCPM_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort),shape=21)+
  scale_color_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#335b76",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)", "HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=13),legend.position = c(0.8,0.9),plot.margin = margin(1,1,1,1, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5)))
dev.off()

######################################### Volcano Plots ######################

######  HC_Conv  -vs- Mild #######
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/All_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Volcano_HC_Conv_vs_Mild.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5,nudge_x = 0.5)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-6.5, 6.5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()


######  HC_Conv  -vs- Severe #####################
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/All_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Volcano_HC_Conv_vs_Severe.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()

######  HC  -vs- Covid #####################
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/All_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Volcano_HC_Conv_vs_COVID.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5,nudge_x = 0.5)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-6.5, 6.5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()

######  Mild  -vs- Severe #####################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/All_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/Volcano_Mild_vs_Severe.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5,nudge_x = 0.5)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()

############################### PIANO KEGG ######################

######  HC_Conv  -vs- Mild #######

library(ComplexHeatmap)
library(circlize)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(10, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/PianoHC_Conv-Mild.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

######  HC_Conv  -vs- Severe #####################
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(15, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/PianoHC_Conv-Severe.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

######  HC  -vs- Covid #####################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(12, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/PianoHC-Covid.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

######  Mild  -vs- Severe #####################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(15, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/NEW_FIGURES/PianoMild-Severe.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()
###############################################

install.packages("/home/anoop/Tools/ComICS_1.0.4.tar.gz", repo=NULL, type="source")

library(ComICS)
data(commons)
data(dcqEx)
results <- dcq(reference_data=immgen_dat, mix_data=lung_time_series_dat, marker_set=DCQ_mar)
names(immgen_dat)

dim(immgen_dat)



library(reshape)
DAT=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Untitled Folder/rna_blood_cell_schmiedel.tsv",sep="\t")
head(DAT)
?cast
Re <- reshape::cast(DAT,fun.aggregate=mean)
dim(Re)
write.table(Re,file="/home/anoop/Desktop/COVID_Omics/EPIC/Untitled Folder/schmiedel_casted.txt",sep="\t",col.names = NA,quote = FALSE)

devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

library(EPIC)
rm(Res)

################################
tdata <- read.delim(file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/TPM_SampleID_GeneName.txt", as.is=T, check.names=F)
bulk <- as.matrix(tdata[,-(1:2)])
rownames(bulk) <- tdata[,1]

rdata <- read.delim(file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Blood_Cells_casted.txt", as.is=T, check.names=F)
REF <- as.matrix(rdata[,-(1:2)])
rownames(REF) <- rdata[,1]
head(REF)

sig=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/NewSigGenes.txt",header = FALSE)
SIG <- as.vector(sig[,1])

refList <- list(refProfiles=REF,sigGenes=SIG)

Res2 <- EPIC(bulk = bulk, reference = refList,scaleExprs=FALSE)

write.table(Res2$mRNAProportions,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/mRNAProp.txt",col.names = NA,quote = FALSE,sep="\t")
write.table(Res2$cellFractions,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/CellFrac.txt",col.names = NA,quote = FALSE,sep="\t")
write.table(Res2$fit.gof,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/OtherRes.txt",col.names = NA,quote = FALSE,sep="\t")

Res3 <- EPIC(bulk = bulk, reference = refList,scaleExprs=FALSE,withOtherCells=FALSE)
write.table(Res3$cellFractions,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/CellFrac2.txt",col.names = NA,quote = FALSE,sep="\t")
?EPIC

library(ggplot2)
library(dplyr)
library(vegan)
CellFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Conv",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/meta.txt",check.names=FALSE,row.names = 1)

permanova <- adonis((CellFr) ~ Group,data = meta, permutations=999, method = "bray")

coef <- coefficients(permanova)["Group1",]
top.coef <- coef[rev(order(abs(coef)))[1:19]]
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Conv.pdf")
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")+theme(plot.margin = margin(1,1,1,3.5, "cm"))
dev.off()
write.table(as.data.frame(sort(top.coef)),file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Conv_Coef.txt",sep="\t",col.names = NA,quote = FALSE)


##############
CellFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Mild.txt",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Meta_HC_Mild.txt",check.names=FALSE,row.names = 1)

permanova <- adonis((CellFr) ~ Group,data = meta, permutations=999, method = "bray")

coef <- coefficients(permanova)["Group1",]
top.coef <- coef[rev(order(abs(coef)))[1:19]]
write.table(as.data.frame(sort(top.coef)),file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Mild_Coef.txt",sep="\t",col.names = NA,quote = FALSE)

#################
CellFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Severe.txt",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Meta_HC_Sev.txt",check.names=FALSE,row.names = 1)

permanova <- adonis((CellFr) ~ Group,data = meta, permutations=999, method = "bray")

coef <- coefficients(permanova)["Group1",]
top.coef <- coef[rev(order(abs(coef)))[1:19]]
write.table(as.data.frame(sort(top.coef)),file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Sev_Coef.txt",sep="\t",col.names = NA,quote = FALSE)

###############
CellFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Mild_Sev.txt",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Met_Mild_Sev.txt",check.names=FALSE,row.names = 1)

permanova <- adonis((CellFr) ~ Group,data = meta, permutations=999, method = "bray")

coef <- coefficients(permanova)["Group1",]
top.coef <- coef[rev(order(abs(coef)))[1:19]]
write.table(as.data.frame(sort(top.coef)),file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Mild_Sev_Coef.txt",sep="\t",col.names = NA,quote = FALSE)

#################


data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MeanFrac.txt")
M=reshape2::melt(data)
head(M)

library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MeanFrac.pdf")
ggplot(M) +
  geom_point(data=M,aes(x=variable,y=factor(CellType,levels = unique(CellType)),size=value,color=variable))+theme_bw()+
  scale_color_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x =element_blank(),
        axis.text.y = element_text(size=10,color="black"),plot.title = element_text(hjust = 0.5),plot.margin = margin(4,4,4.5,3.5, "cm"),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),
        legend.direction = "vertical",
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Mean Cell Proportion",ncol = 1),color=guide_legend(title = NULL))
dev.off()

############################
library(ComplexHeatmap)
library(circlize)
CelFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/CellFrac.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(CelFr)
col_fun2 = colorRamp2(c(1,0.3,0.25,0.2,0.15,0.1,0.05,0), c("#2a4e6c","#386890","#386890","#3f75a2","#4682b4","#588ebb","#7da7ca","#ecf2f7"))
sample=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Meta.txt",row.names = 1)

ha = HeatmapAnnotation(df = sample,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10), 
                                                                   labels_gp = gpar(fontsize = 10))),
                       col = list(Group=c("HC"="#48690E","Mild"="#ffd700","Severe"="#ffa500")))

H1=Heatmap(as.matrix(t(CelFr)),cluster_columns = TRUE,show_column_names = FALSE,col = col_fun2,top_annotation = ha,
           column_dend_height = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),column_split =sample$Group,cluster_rows = FALSE,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Cell proportion",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(7, "cm"),width  = unit(8, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/CellFrac.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

cof=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Mild_Coef.txt")
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Mild_Coef.pdf")
ggplot(data=cof, aes(x=Coef, y=factor(CellType,levels = unique(CellType)))) +labs(x="Permanova coefficient",title = "HC-vs-Mild (pval=0.001)")+
  geom_bar(stat="identity",fill="steelblue")+theme(plot.margin = margin(4.5,4.5,4.5,3.5, "cm"),plot.title = element_text(size=10,colour = "black"),
                                                   axis.title.y =element_blank(),axis.title.x = element_text(size=10,colour = "black"),
                                  axis.text = element_text(size=10,colour = "black"))
dev.off()

cof=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Sev_Coef.txt")
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_Sev_Coef.pdf")
ggplot(data=cof, aes(x=Coef, y=factor(CellType,levels = unique(CellType)))) +labs(x="Permanova coefficient",title = "HC-vs-Severe (pval=0.001)")+
  geom_bar(stat="identity",fill="#817a99")+theme(plot.margin = margin(4.5,4.5,4.5,3.5, "cm"),axis.title.y =element_blank(),
                                                 axis.title.x = element_text(size=10,colour = "black"),plot.title = element_text(size=10,colour = "black"),
                                                   axis.text = element_text(size=10,colour = "black"))
dev.off()


cof=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Mild_Sev_Coef.txt")
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/Mild_Sev_Coef.pdf")
ggplot(data=cof, aes(x=Coef, y=factor(CellType,levels = unique(CellType)))) +labs(x="Permanova coefficient",title = "Mild-vs-Severe (pval=0.058)")+
  geom_bar(stat="identity",fill="#41929b")+theme(plot.margin = margin(4.5,4.5,4.5,3.5, "cm"),axis.title.y =element_blank(),plot.title = element_text(size=10,colour = "black"),
                                                 axis.title.x = element_text(size=10,colour = "black"),
                                                 axis.text = element_text(size=10,colour = "black"))
dev.off()


##########

CellFr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/CellFrac.txt",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Meta2.txt",check.names=FALSE,row.names = 1)

permanova <- adonis((CellFr) ~ Group,data = meta, permutations=999, method = "bray")

coef <- coefficients(permanova)["Group1",]
top.coef <- coef[rev(order(abs(coef)))[1:19]]
write.table(as.data.frame(sort(top.coef)),file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_COVID_Coef.txt",sep="\t",col.names = NA,quote = FALSE)


cof=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_COVID_Coef.txt")
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/HC_COVID_Coef.pdf")
ggplot(data=cof, aes(x=Coef, y=factor(CellType,levels = unique(CellType)))) +labs(x="Permanova coefficient",title = "HC-vs-Covid (pval=0.001)")+
  geom_bar(stat="identity",fill="#524cff")+theme(plot.margin = margin(4.5,4.5,4.5,3.5, "cm"),axis.title.y =element_blank(),plot.title = element_text(size=10,colour = "black"),
                                                 axis.title.x = element_text(size=10,colour = "black"),
                                                 axis.text = element_text(size=10,colour = "black"))
dev.off()

###########################

data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MedianCellFr.txt")
M=reshape2::melt(data)
head(M)
#scale_color_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MedianCellFr.txt.pdf")
ggplot(M) +
  geom_point(data=M,aes(x=variable,y=factor(CellType,levels = unique(CellType)),size=value,color=value))+theme_bw()+
  scale_color_gradient(low="#669cc0",high="#00365a")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x =element_text(size=10,color="black",angle =45, vjust =1.2, hjust=1),
        axis.text.y = element_text(size=10,color="black"),plot.title = element_text(hjust = 0.5),plot.margin = margin(5,4,5,3.5, "cm"),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),
        legend.direction = "vertical",
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Median Cell Proportion",ncol = 1),color=guide_legend(title = NULL))+
  scale_size_continuous(range = c(2,5))
dev.off()


pval=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/Kruskal_Wallis_Result.tab",header = T,check.names = FALSE)
head(pval)
pval$AdjPval_KW =p.adjust(pval$KW,method = "BH")
pval$AdjPval_Mild_HC =p.adjust(pval$`Mild-HC`,method = "BH")
pval$AdjPval_Sev_HC =p.adjust(pval$`Severe-HC`,method = "BH")
pval$AdjPval_Mild_Sev =p.adjust(pval$`Severe-Mild`,method = "BH")
write.table(pval,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/Kruskal_Wallis_Result_Adj.tab",sep="\t",col.names = NA, quote = FALSE)

pval=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/pair.tab",header = T,check.names = FALSE)
head(pval)
pval$AdjPval =p.adjust(pval$pval,method = "BH")

write.table(pval,file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/pairAdj.txt",sep="\t",col.names = NA, quote = FALSE)


box=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/CellFr.tab",header = TRUE,check.names = FALSE)

my_comp <- list( c("HC", "Mild"), c("HC", "Severe"), c("Mild", "Severe") )

P7=ggplot(box,aes(x=Group,y=classical_monocyte,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Classical monocytes")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)
?stat_compare_means
P1=ggplot(box,aes(x=Group,y=neutrophil,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Neutrophil")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P3=ggplot(box,aes(x=Group,y=`MAIT_T-cell`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "MAIT T-cell")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P4=ggplot(box,aes(x=Group,y=`gdT-cell`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "gd T-cell")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P6=ggplot(box,aes(x=Group,y=`intermediate_monocyte_2`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Intermediate monocyte")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P2=ggplot(box,aes(x=Group,y=`memory_B-cell`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Memory B-cell")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P5=ggplot(box,aes(x=Group,y=`memory_CD8_T-cell`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","HC (CoV-2 Ab+)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Memory CD8 T-cell")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

P8=ggplot(box,aes(x=Group,y=`non-classical_monocyte`,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                     values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(labels = c("Healthy Control (HC)","Hospitalized-mild","Hospitalized-severe"),
                    values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+labs(title = "Non-classical monocyte")+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.text = element_text(size = 9),axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),legend.title = element_blank(),
        panel.border = element_blank(),legend.position = c(1.7,0.5))+
  stat_compare_means(comparisons = my_comp,hide.ns=TRUE,size=0)

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/SigniCells.pdf",height = 5,width =10)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,nrow = 2,ncol = 5)+theme(plot.margin = margin(2,2,2,0.5, "cm"))
dev.off()


############


DAT=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/Comm.txt")

Re <- reshape::cast(DAT,fun.aggregate=mean)
head(Re)

write.table(Re,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/AllCorr_Casted.txt",sep = "\t",col.names = NA,quote = FALSE)

com=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/Top5.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Comm=c("CC1_1"="#6666FF","CC1_2"="#FFCC33","CC1_3"="#669900",
                                     "CC1_4"="#CC6600","CC1_5"="#00CCFF","CC1_6"="#00CCCC","CC1_7"="#999900","CC1_8"="#FF6666","CC1_9"="#FFCC66","CC1_10"="#0000FF")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                        grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Comm=c("CC1_1"="#6666FF","CC1_2"="#FFCC33","CC1_3"="#669900",
                                          "CC1_4"="#CC6600","CC1_5"="#00CCFF","CC1_6"="#00CCCC","CC1_7"="#999900","CC1_8"="#FF6666","CC1_9"="#FFCC66","CC1_10"="#0000FF")))


col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))


corr=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/TopCorr_Casted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha3,left_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(7, "cm"),width  = unit(7, "cm"))



pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/TopCorr_Casted.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


##############
dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C2/Ip.txt",header = TRUE)
head(dat)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C2/C2.pdf")
ggplot(data) +
  geom_point(data=dat,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#ffe084",high="#997a1e")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(5,4.5,7,2, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()


dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C6/Ip.txt",header = TRUE)
head(dat)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C6/C6.pdf")
ggplot(data) +
  geom_point(data=dat,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#66e0e0",high="#006666")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(6,4.5,7,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()


dat=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C7/Ip.txt",header = TRUE)
head(dat)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/NEW/results/KEGG_Meta/C7/C7.pdf")
ggplot(data) +
  geom_point(data=dat,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#c1c166",high="#5b5b00")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(7,3,7,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

############


library(piano)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/All_Result_2.txt",header = TRUE,row.names = 1)
p <- data[7]
lfc<-data[4]
head(p)
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/HumanGEM_subsystem_GSC.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,geneSetStat="reporter",
                 nPerm = 1000)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/ReporterSubSystem.txt",sep="\t",col.names = NA)
?runGSA

######################### MUBR SNF Clusters

library(doParallel)
library(MUVR)

nCore=10
nRep=25
nOuter=8
varRatio=0.8
method='RF'
data=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/Mild_Severe/XX.txt",check.names = FALSE)
data$METs <- NULL

write.table(data,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/Mild_Severe/XXX.txt",sep="\t",row.names = FALSE,quote = FALSE)

XX=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/test.txt",check.names = FALSE)
ncol(XX)
Y=c(rep("Mild",26),rep("Severe",11))
YY=as.factor(Y)

cl=makeCluster(nCore)
registerDoParallel(cl)
classModel = MUVR(X=XX, Y=YY, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)
cbind(YY, classModel$yClass)
classModel$miss
classModel$nVar

vip=getVIP(classModel, model='min')

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/Model_Validation.pdf")
plotVAL(classModel)
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/SwimLanePlot.pdf")
plotMV(classModel, model='min')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/StabilityPlot.pdf")
plotStability(classModel, model='mid')
dev.off()

pdf("/home/anoop/Desktop/COVID_Omics/METABOLOME/NEW_Plan/MUVR/VIP.pdf")
plotVIP(classModel, model='mim')
dev.off()

vip
write.table(vip,file="/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/Variables.txt",sep="\t",col.names = NA,quote = FALSE)

######################################### BioMarkers SNF Clusters

Bean=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/BIOMARK.txt",check.names = FALSE)

p1=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`eicosanedioate (C20-DC)`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "Eicosanedioate (C20-DC)")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))

p2=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`mannose`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "Mannose")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))

p3=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`eicosenedioate (C20:1-DC)*`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "Eicosenedioate (C20:1-DC)*")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))

p4=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`hydantoin-5-propionate`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "Hydantoin-5-propionate")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))

p5=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`6-bromotryptophan`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-bromotryptophan")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))

p6=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`1-palmitoyl-2-oleoyl-GPC (16:0||18:1)`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "1-palmitoyl-2-oleoyl-GPC (16:0||18:1)")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.7,size = 9))


p7=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`4-hydroxyphenylacetate`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "4-hydroxyphenylacetate")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))


p8=ggplot(Bean,aes(x=factor(METs,levels = unique(METs)),y=log2(`6-oxopiperidine-2-carboxylate`),fill=METs,color=METs))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(measurement)",title = "6-oxopiperidine-2-carboxylate")+
  theme_bw()+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400",SNF_Cls4_HC="#3f75a2",SNF_Cls1_Mild="#8b814c",SNF_Cls3_Severe="#cd5555"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500",SNF_Cls4_HC="#4682b4",SNF_Cls1_Mild="#cdbe70",SNF_Cls3_Severe="#ee6363"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(-0.3,-0.1),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 8),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 9))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1)) 

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/MUVR/BIOMARKers.pdf",height = 5,width =10)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow = 2,ncol = 4)+theme(plot.margin = margin(0.5,0.5,1,0.5, "cm"))
dev.off()



pval=read.delim("/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/CoExpression/Proteom/results/Species/CommCheck/CC1_4/CUSTOM140713934812256.ART.enrichr.reports.txt",header = TRUE)
head(pval)
pval$AdjPval =p.adjust(pval$P.value,method = "BH")
write.table(pval,file="/home/anoop/Desktop/Proteo-Transcriptomics/New_Analysis/PLS/PLS_NEW/CoExpression/Proteom/results/Species/CommCheck/CC1_4/ResultsAdjested.tab",sep="\t",col.names = NA, quote = FALSE)


count=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Combined/GSEA/Coding_Count_Hg38.txt",header = TRUE,row.names = 1)

library(NOISeq)
tmm=tmm(as.matrix(count))
head(tmm)

write.table(tmm,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Combined/GSEA/Tmm.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
data=read.delim("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Combined/results/Species/CommCheck/CC3_3_AllKegg/Ip.txt",header = TRUE)
head(data)
pdf("/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Combined/results/Species/CommCheck/CC3_3_AllKegg/C3.pdf")
ggplot(data) +
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#84b987",high="#0a740f")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(3.5,5,3.5,2, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

#################

################################# Pinao


library(piano)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/All_Result_2.txt",header = TRUE,row.names = 2)
p <- data[7]
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/PianoKEGG.txt",sep="\t",col.names = NA)


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/All_Result_2.txt",header = TRUE,row.names = 2)
p <- data[7]
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/PianoKEGG.txt",sep="\t",col.names = NA)


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/All_Result_2.txt",header = TRUE,row.names = 2)
p <- data[7]
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/PianoKEGG.txt",sep="\t",col.names = NA)


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/All_Result_2.txt",header = TRUE,row.names = 2)
p <- data[7]
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/PianoKEGG.txt",sep="\t",col.names = NA)

######################################### Piano

library(ComplexHeatmap)
library(circlize)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a)p adj (mix.dir.up)"="#a6a6a6",
                                          "b)p adj (mix.dir.dn)"="#8c8c8c","c)p adj (non-dir.)"="#808080",
                                          "d)p adj (dist.dir.dn)"="#666666","e)p adj (dist.dir.up)"="#4c4c4c")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(10, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/PianoKEGG.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

######################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a)p adj (mix.dir.up)"="#a6a6a6",
                                          "b)p adj (mix.dir.dn)"="#8c8c8c","c)p adj (non-dir.)"="#808080",
                                          "d)p adj (dist.dir.dn)"="#666666","e)p adj (dist.dir.up)"="#4c4c4c")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(15,"cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Severe/PianoKEGG.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

######################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a)p adj (mix.dir.up)"="#a6a6a6",
                                          "b)p adj (mix.dir.dn)"="#8c8c8c","c)p adj (non-dir.)"="#808080",
                                          "d)p adj (dist.dir.dn)"="#666666","e)p adj (dist.dir.up)"="#4c4c4c")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(15,"cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Covid/PianoKEGG.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

##########################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/pianokeggIP.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a)p adj (mix.dir.up)"="#a6a6a6",
                                          "b)p adj (mix.dir.dn)"="#8c8c8c","c)p adj (non-dir.)"="#808080",
                                          "d)p adj (dist.dir.dn)"="#666666","e)p adj (dist.dir.up)"="#4c4c4c")))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(15,"cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Mild-Severe/PianoKEGG.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

#######################################

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res/HMR.txt",row.names = 1,header = TRUE,check.names = FALSE)
col_fun1 = colorRamp2(c(-1,0,1), c("#005900" ,"white","#b20000"))

H1=Heatmap(as.matrix(DAT),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_heatmap_legend = FALSE,
           left_annotation = rowAnnotation(Subsytem = pwy$Subsystem,col=colours,show_legend=TRUE,show_annotation_name=FALSE),row_title_gp = gpar(fontsize=10),
           show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(20,"cm"),width  = unit(4, "cm"))

pwy=read.delim("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res/subs.txt",row.names = 1)
head(pwy)
ha = rowAnnotation(df = pwy,show_annotation_name = TRUE,)

pdf("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res/HMR.pdf",width = 7,height=10)
draw(H1, merge_legend = TRUE)
dev.off()

########

DAT=read.delim("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res2/fba.txt",row.names = 1,header = TRUE,check.names = FALSE)
col_fun1 = colorRamp2(c(-1,0,1), c("#005900" ,"white","#b20000"))

H1=Heatmap(as.matrix(DAT),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_heatmap_legend = FALSE,
           left_annotation =ha,
           show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(15,"cm"),width  = unit(4, "cm"))

pwy=read.delim("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res2/Subs.txt",row.names = 1)
head(pwy)
ha = rowAnnotation(df = pwy,show_annotation_name = TRUE)

pdf("/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/Analyis/res2/HMR.pdf",width = 7,height=10)
draw(H1, merge_legend = TRUE)
dev.off()

###########


DAT=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Analysis/res/fbaData.txt",row.names = 1,header = TRUE,check.names = FALSE)
col_fun1 = colorRamp2(c(-1,0,1), c("#005900" ,"white","#b20000"))

H1=Heatmap(as.matrix(DAT),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_heatmap_legend = FALSE,
           left_annotation =ha,
           show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(18,"cm"),width  = unit(4, "cm"))

pwy=read.delim("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Analysis/res/sbs.txt",row.names = 1)
head(pwy)
ha = rowAnnotation(df = pwy,show_annotation_name = TRUE)

pdf("/home/anoop/Desktop/COVID_Omics/COVID_network_fusion/Analysis/res/HMR.pdf",width = 7,height=10)
draw(H1, merge_legend = TRUE)
dev.off()

###########

Data=read.delim("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/DC.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list( c("HC", "Mild"),c("HC", "Severe"),c("Mild", "Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/DC.pdf",width = 17.5,height = 8)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 8,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.1,-0.1),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(5,0.1,5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/MDSc.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list( c("HC", "Mild"),c("HC", "Severe"),c("Mild", "Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/MDSC.pdf",width = 17.5,height = 8)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.1,-0.1),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(2,0.1,2,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()



Data=read.delim("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/Mono.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list( c("HC", "Mild"),c("HC", "Severe"),c("Mild", "Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/Cells/MW/Plot/Monocyte.pdf",width = 17.5,height = 10)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 6,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.1,-0.05),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(2,0.1,2,0.1, "cm"),strip.text.x = element_text(size = 11.5),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()

rm(list=ls())
packageVersion("piano")
#########
install.packages("mlbench")
library(tidyverse)
library(caret)
data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)
sample_n(PimaIndiansDiabetes2, 3)
# Split the data into training and test set
set.seed(123)
training.samples <- PimaIndiansDiabetes2$diabetes %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- PimaIndiansDiabetes2[training.samples, ]
test.data <- PimaIndiansDiabetes2[-training.samples, ]

model <- glm( diabetes ~., data = train.data, family = binomial)
summary(model)
probabilities <- model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
mean(predicted.classes == test.data$diabetes)
predicted.classes

model <- glm( diabetes ~ mass, data = train.data, family = binomial)
summary(model)$coef


dat=read.delim("/home/anoop/Desktop/COVID_Omics/Logistic_regression/Facs.txt",row.names = 1)
Lof=log2(dat)
head(Lof)
write.table(Lof,file="/home/anoop/Desktop/COVID_Omics/Logistic_regression/FacsLog2.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Logistic_regression/Olink.txt",row.names = 1)
Lof=log2(dat)
head(dat)
write.table(Lof,file="/home/anoop/Desktop/COVID_Omics/Logistic_regression/OlinkLog2.txt",sep="\t",col.names = NA,quote = FALSE)

#######

Data=read.delim("/home/anoop/Desktop/COVID_Omics/Logistic_regression/Olink.txt",row.names = 1,check.names = FALSE)
set.seed(123)
training.samples <- Data$Covid %>% createDataPartition(p = 0.8, list = FALSE)
head(Data)

train.data  <- Data[training.samples, ]
test.data <- Data[-training.samples, ]

head(train.data)
dim(test.data)

model <- glm( Covid ~., data = train.data, family = binomial)
summary(model)
?summary
(res$coefficients)
head(as.data.frame(model$effects))
coef(summary(model))[,'Pr(>|z|)']
XX=cor(Data[-1])
head(XX)
library(reshape2)
MM=melt(XX)
head(MM)
write.table(MM,file="/home/anoop/Desktop/COVID_Omics/Logistic_regression/test.txt",sep="\t",col.names = NA,quote = FALSE)
DD=vif(model)
library(car)


############

install.packages("Boruta")

data <- read.delim('/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Input_FACS.tab',header=TRUE,sep='\t')
head(data)
library(Boruta)

set.seed(123)
boruta_res <- Boruta(Covid~., data = data, doTrace = 2,mtry=20,ntree=1000,maxRuns=500)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/FACS_importance.pdf",width = 10,height = 5)
plot(boruta_res, las = 2, cex.axis = 0.5)
dev.off()
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/FACS_HistoryImp.pdf",width = 10,height = 5)
plotImpHistory(boruta_res)
dev.off()
library(caret)
set.seed(123)
ind <- sample(2, nrow(data), replace = T, prob = c(0.6, 0.4))
train <- data[ind==1,]
head(train)
test <- data[ind==2,]
library(randomForest)
set.seed(456)
rf_all <- randomForest(Covid~., data = train)
p <- predict(rf_all, test)
confusionMatrix(p, test$Covid)

getNonRejectedFormula(boruta_res)
set.seed(789)

rf_boruta <- randomForest(Covid ~ cDC.1 + cDC.2.CCR5 + pDC + pDC.CCR2 + pDC.CCR5 + Monocyte.derived.CD141.CCR2 + 
                            Monocyte.derived.CD141.CCR5 + Monocyte.derived.CD1c + G.MDSC + 
                            G.MDSC.CCR2 + G.MDSC.CX3CR1 + M.MDSC + M.MDSC.CCR5 + CM + 
                            IM.CX3CR1 + NCM + NCM.CCR5.CCR1 + NCM.CCR5, data=train)

p <- predict(rf_boruta, test)
confusionMatrix(p, test$Covid)
boruta_table <- attStats(boruta_res)
write.table(boruta_table,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/FACS_features.txt",sep="\t",col.names = NA,quote = FALSE)

rf_prediction <- predict(rf_boruta,test,type = "prob")
ROC_rf <- roc(test$Covid, rf_prediction[,2])
ROC_rf_auc <- auc(ROC_rf)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/FACS_ROC.pdf")
plot(ROC_rf)
dev.off()

##################

data <- read.delim('/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Metabolites.txt',header=TRUE,sep='\t',row.names = 1)
head(data)
library(Boruta)

set.seed(123)
boruta_res <- Boruta(Covid~., data = data, doTrace = 2,mtry=20,ntree=1000,maxRuns=500)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Mets_importance.pdf",width = 10,height = 7)
plot(boruta_res, las = 2, cex.axis = 0.5)
dev.off()
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Mets_HistoryImp.pdf",width = 10,height = 5)
plotImpHistory(boruta_res)
dev.off()

set.seed(123)
ind <- sample(2, nrow(data), replace = T, prob = c(0.6, 0.4))
train <- data[ind==1,]
head(train)
test <- data[ind==2,]
library(randomForest)
set.seed(456)
rf_all <- randomForest(Covid~., data = train)
p <- predict(rf_all, test)
confusionMatrix(p, test$Covid)

getNonRejectedFormula(boruta_res)
set.seed(789)

rf_boruta <- randomForest(Covid ~ X.N.1....N.8...acetylspermidine_Amino.Acid + X1.ribosyl.imidazoleacetate._Amino.Acid + 
                            X2.ketocaprylate_Amino.Acid + X2.oxoarginine._Amino.Acid + 
                            X4.hydroxyglutamate_Amino.Acid + X5.hydroxyindole.sulfate_Amino.Acid + 
                            X6.bromotryptophan_Amino.Acid + alanine_Amino.Acid + catechol.glucuronide_Amino.Acid + 
                            cysteine.s.sulfate_Amino.Acid + dopamine.4.sulfate_Amino.Acid + 
                            homovanillate..HVA._Amino.Acid + imidazole.propionate_Amino.Acid + 
                            indolelactate_Amino.Acid + indolepropionate_Amino.Acid + 
                            methylsuccinate_Amino.Acid + methylsuccinoylcarnitine_Amino.Acid + 
                            N.acetyltryptophan_Amino.Acid + N2.N5.diacetylornithine_Amino.Acid + 
                            oxindolylalanine_Amino.Acid + S.carboxyethylcysteine_Amino.Acid + 
                            hydroxyproline_Amino.Acid + tryptophan_Amino.Acid + tyramine.O.sulfate_Amino.Acid + 
                            tyrosine_Amino.Acid + mannose_Carbohydrate + N.acetylneuraminate_Carbohydrate + 
                            bilirubin..E.Z.or.Z.E.._Cofactors.and.Vitamins + carotene.diol..1._Cofactors.and.Vitamins + 
                            retinol..Vitamin.A._Cofactors.and.Vitamins + citraconate..glutaconate_Energy + 
                            citrate_Energy + X1.2.dilinoleoyl.GPC..18.2..18.2._Lipid + 
                            X1.lignoceroyl.GPC..24.0._Lipid + X1.myristoyl.2.arachidonoyl.GPC..14.0..20.4.._Lipid + 
                            X1.palmitoleoyl.2.linolenoyl.GPC..16.1..18.3.._Lipid + X11beta.hydroxyandrosterone.glucuronide_Lipid + 
                            X2.aminoheptanoate_Lipid + X2.aminooctanoate_Lipid + X2.hydroxydecanoate_Lipid + 
                            X2.hydroxyoctanoate_Lipid + X2.hydroxysebacate_Lipid + X3.carboxy.4.methyl.5.pentyl.2.furanpropionate..3.CMPFP._Lipid + 
                            andro.steroid.monosulfate.C19H28O6S..1.._Lipid + androstenediol..3beta.17beta..disulfate..2._Lipid + 
                            androsterone.glucuronide_Lipid + arachidonoyl.ethanolamide_Lipid + 
                            docosadioate..C22.DC._Lipid + eicosanedioate..C20.DC._Lipid + 
                            glycerol.3.phosphate_Lipid + glycosyl.ceramide..d18.1..20.0..d16.1..22.0.._Lipid + 
                            N.acetyl.2.aminooctanoate._Lipid + N.stearoyl.sphingosine..d18.1..18.0.._Lipid + 
                            oleoyl.ethanolamide_Lipid + picolinoylglycine_Lipid + sphingomyelin..d17.1..14.0..d16.1..15.0.._Lipid + 
                            tetrahydrocortisone.glucuronide..5._Lipid + X5.methyluridine..ribothymidine._Nucleotide + 
                            adenosine.3..5..cyclic.monophosphate..cAMP._Nucleotide + 
                            allantoin_Nucleotide + gamma.glutamylalanine_Peptide + gamma.glutamylcitrulline._Peptide + 
                            gamma.glutamylhistidine_Peptide + gamma.glutamylleucine_Peptide + 
                            gamma.glutamyltyrosine_Peptide + prolylglycine_Peptide + 
                            valylglycine_Peptide + valylleucine_Peptide, data=train)

p <- predict(rf_boruta, test)
confusionMatrix(p, test$Covid)
boruta_table <- attStats(boruta_res)
write.table(boruta_table,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Mets_features.txt",sep="\t",col.names = NA,quote = FALSE)

rf_prediction <- predict(rf_boruta,test,type = "prob")
ROC_rf <- roc(test$Covid, rf_prediction[,2])
ROC_rf$sensitivities
ROC_rf$specificities
ROC_rf_auc <- auc(ROC_rf)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Mets_ROC.pdf")
plot(ROC_rf)
dev.off()

###################################

data <- read.delim('/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Olink.txt',header=TRUE,sep='\t',row.names = 1)
head(data)
library(Boruta)

set.seed(123)
boruta_res <- Boruta(Covid~., data = data, doTrace = 2,mtry=20,ntree=1000,maxRuns=500)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Olink_importance.pdf",width = 10,height = 7)
plot(boruta_res, las = 2, cex.axis = 0.5)
dev.off()
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Olink_HistoryImp.pdf",width = 10,height = 5)
plotImpHistory(boruta_res)
dev.off()

set.seed(123)
ind <- sample(2, nrow(data), replace = T, prob = c(0.6, 0.4))
train <- data[ind==1,]
head(train)
test <- data[ind==2,]
library(randomForest)
detach("package:randomForest", unload=TRUE)
set.seed(456)
rf_all <- randomForest(Covid~., data = train)
p <- predict(rf_all, test)
confusionMatrix(p, test$Covid)

getNonRejectedFormula(boruta_res)
set.seed(789)

rf_boruta <- randomForest(Covid ~ IL8 + MCP.3 + IL6 + CXCL11 + Gal.9 + CD40 + TNFSF14 + 
                            CD28 + CCL19 + MCP.2 + CCL4 + IL15 + PD.L1 + HGF + CXCL10 + 
                            IL10 + CCL3 + PTN + CXCL12 + IFN.gamma + CASP.8 + CXCL13 + 
                            VEGFA + IL12RB1 + TNF + GZMB + CSF.1, data=train)

p <- predict(rf_boruta, test)
confusionMatrix(p, test$Covid)
boruta_table <- attStats(boruta_res)
write.table(boruta_table,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Olink_features.txt",sep="\t",col.names = NA,quote = FALSE)

rf_prediction <- predict(rf_boruta,test,type = "prob")
ROC_rf <- roc(test$Covid, rf_prediction[,2])
ROC_rf$sensitivities
ROC_rf$specificities
ROC_rf_auc <- auc(ROC_rf)
pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Olink_ROC.pdf")
plot(ROC_rf)
dev.off()
############


Data=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/cDC-1 analysis.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/ubha/cDC-1 analysis.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("HIV"="#964C44","Ctrl"="#726F82"))+
  scale_fill_manual(values=c("HIV"="#FB8072","Ctrl"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()

Data=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/cDC-2 analysis.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/ubha/cDC-2 analysis.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("HIV"="#964C44","Ctrl"="#726F82"))+
  scale_fill_manual(values=c("HIV"="#FB8072","Ctrl"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/iDC analysis.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("//home/anoop/Desktop/COVID_Omics/ubha/iDC analysis.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("HIV"="#964C44","Ctrl"="#726F82"))+
  scale_fill_manual(values=c("HIV"="#FB8072","Ctrl"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/COVID_Omics/ubha/pDC analysis.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/ubha/pDC analysis.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("HIV"="#964C44","Ctrl"="#726F82"))+
  scale_fill_manual(values=c("HIV"="#FB8072","Ctrl"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


#################

Data=read.delim("/home/anoop/Desktop/COVID_Omics/Cells/NEW_DATA/s7.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list( c("HC", "Mild"),c("HC", "Severe"),c("Mild", "Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/Cells/NEW_DATA/set7.pdf",width = 17.5,height = 10)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 5,scales ="free" )+
  geom_jitter(shape=1, size=2.5,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.1,-0.05),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,4,4.5,4, "cm"),strip.text.x = element_text(size = 11.5),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()
packageVersion("ComplexHeatmap")

####################

library(psych)
Genes=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Correlation/Merged.txt",check.names=FALSE,header = T)
Res=corr.test(as.matrix(Genes),use = "pairwise",method="spearman",adjust="BH")

head(Res$r)
corr=melt(Res$r)
head(corr)
pval=melt(Res$p)

corr$G1=pval$Var1
corr$G2=pval$Var2
corr$pvalue=pval$value

write.table(corr,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/Correlation/CorrResults.txt",sep="\t",col.names = NA,quote = FALSE)



######################

Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC1.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC1.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC2.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC2.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC3.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/DC3.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/MDSC.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/MDSC.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/Mono1.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/Mono1.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


Data=read.delim("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/Mono2.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/Shubha/CoCoMo_FACS_data/Mono2.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 7,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("ART"="#964C44","HC"="#726F82"))+
  scale_fill_manual(values=c("ART"="#FB8072","HC"="#BEBADA"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(4.5,0.1,4.5,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


#####################

library(ComplexHeatmap)
library(circlize)

snf=read.delim("/home/anoop/Desktop/COVID_Omics/Personal_GEM/Essentiality/Analysis/Cls_Info.txt",row.names = 1)
head(snf)

ha = HeatmapAnnotation(df = snf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"), height = unit(3, "cm"),
                       annotation_legend_param = list(Cluster = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                  grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize =8))),
                       col = list(Cluster=c("cls1"="#8DD3C7","cls2"="#FB8072","cls3"="#BEBADA","cls4"="#666666")))

clini=read.delim("/home/anoop/Desktop/COVID_Omics/Personal_GEM/Essentiality/Analysis/Clini_Info.txt",row.names = 1)
head(clini)

ha2 = HeatmapAnnotation(df = clini,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"), height = unit(3, "cm"),
                       annotation_legend_param = list(Cohort = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                  grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8))),
                       col = list(Cohort=c("Mild"="#ffd700","Severe"="#ffa500","Healthy_Conv"="#8BB443","Healthy"="#48690E")))

data=read.delim("/home/anoop/Desktop/COVID_Omics/Personal_GEM/Essentiality/Analysis/Covid_Samples.txt",row.names = 1)
col_fun1 = colorRamp2(c(1,0), c("#3e5a70","white"))

H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(0.5, "mm"),show_heatmap_legend = FALSE,
           top_annotation  =ha,bottom_annotation = ha2,heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8),
                                                                                  labels_gp = gpar(fontsize = 8)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 4),height  = unit(12, "cm"),width  = unit(10, "cm"),
           column_split =snf$Cluster)
pdf("/home/anoop/Desktop/COVID_Omics/Personal_GEM/Essentiality/Analysis/Covid_Samples.pdf")
draw(H1,merge_legend = TRUE)
dev.off()
?tiff
tiff("/home/anoop/Desktop/COVID_Omics/Personal_GEM/Essentiality/Analysis/Covid_Samples.tiff",width = 7, height = 7, units = 'in', res = 500)
draw(H1,merge_legend = TRUE)
dev.off()

?Heatmap
#########

IL=read.delim("/home/anoop/Desktop/COVID_Omics/Olink/IL_data.txt",row.names = 1)
head(IL)

XX=scale(IL)
head(XX)

write.table(XX,file="/home/anoop/Desktop/COVID_Omics/Olink/IL_data_scaled.txt",sep="\t",col.names = NA,quote = FALSE)

data=read.delim("/home/anoop/Desktop/COVID_Omics/Olink/Median.txt")
library(ggpubr)

pdf("/home/anoop/Desktop/COVID_Omics/Olink/Median.pdf")
ggplot(data,aes(x=factor(Group,levels = unique(Group)),y=IL,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("Cluster1"="#547e77","Cluster2"="#af594f","Cluster3"="#726f82","Cluster4"="#3d3d3d"))+
  scale_fill_manual(values=c("Cluster1"="#8DD3C7","Cluster2"="#FB8072","Cluster3"="#BEBADA","Cluster4"="#666666"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "right",
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 10),
        plot.margin = margin(4.5,3.5,5,3.5, "cm"),strip.text.x = element_text(size = 10),
        axis.text.y = element_text(color="black",size=10),plot.title = element_text(hjust = 0.5,size = 10))+
  guides(color = guide_legend(ncol = 1),fill = guide_legend(ncol = 1))
dev.off()

########

library(psych)
library(reshape2)
cd=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/HC_cd.txt",check.names=FALSE,header = TRUE)
facs=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/HC_facs.txt",check.names=FALSE,header = TRUE)
Corr=corr.test(as.matrix(cd),as.matrix(facs), use = "pairwise",method="spearman",adjust="BH")
write.table(Corr$r,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/HC_R.txt",sep="\t",col.names = NA,quote = FALSE)


Corr2=corr.test(as.matrix(cd),as.matrix(facs), use = "pairwise",method="spearman",adjust="none")

head(Corr$r)
corr=melt(Corr$r)
head(corr)
pval=melt(Corr$p)
rowP=melt(Corr2$p)
head(rowP)
corr$BH=pval
corr$Nominal=rowP
write.table(corr,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/Mild_Results.txt",sep="\t",col.names = NA,quote = FALSE)

#######
library(ComplexHeatmap)
library(circlize)
info=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/Info.txt",row.names = 1)
head(info)
ha = HeatmapAnnotation(df = info,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Type=list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8),
                                                                 Cohort=list(direction = "horizontal",grid_width = unit(0.5, "cm"))),
                                                                           grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                        col = list(Type=c("CD4"="#ccae94","CD8"="#86959f"),Cohort=c("Mild"="#ffd700","Severe"="#ffa500","HC"="#48690E")))


col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))


corr=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/Correlation.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),column_split = info$Cohort,row_names_side = "left",
           name = "Correlation Coeficient",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(8, "cm"),width  = unit(5, "cm"))



pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/Correlation.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/TCA_Log2.txt")
head(data)
library(ggpubr)
ggbarplot(ToothGrowth, x = "Group", y = "malate", 
          add = c("mean_se", "jitter"),
          color = "supp", palette = c("#00AFBB", "#E7B800"),
          position = position_dodge(0.8))


#############


Data=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/Ip_Plot.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data,id.vars = "Group")
head(M)
library(ggpubr)
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/Ip_Plot.pdf",width = 19,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 3,scales ="free" )+
  geom_jitter(shape=1, size=5,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c("COVID"="#cc8400","HC"="#39540b"))+
  scale_fill_manual(values=c("COVID"="#ffa500","HC"="#48690E"))+
  theme_bw()+stat_compare_means(label = "p.format")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(0.05,-0.03),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(7,8,7,8, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))
dev.off()


dat=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP.txt",row.names = 1,check.names = FALSE)
head(dat)
XX=dat*31.73
head(XX)

write.table(XX,file="/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP_Multiplied.txt",sep="\t",col.names = NA,quote = FALSE)

#####################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP_Multiplied.txt",header = TRUE,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP_Multiplied.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP_Multiplied_Z.txt",sep="\t",quote = FALSE,col.names = NA)


clini=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/Info.txt",row.names = 1)
head(clini)

ha2 = HeatmapAnnotation(df = clini,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"), height = unit(3, "cm"),
                        annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                     grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8))),
                        col = list(Group=c("Mild"="#ffd700","Severe"="#ffa500","HC"="#48690E")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#5D5DFF","#4C4CFF","#0000CC"))
data=read.delim("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/IP_Multiplied_Z.txt",row.names = 1)

head(data)
H1=Heatmap(as.matrix((data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(0.5, "mm"),
           top_annotation  =ha2,heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8),
                                                                                  labels_gp = gpar(fontsize = 8)),
           name = "Z-Scaled_Conc",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(7, "cm"),width  = unit(10, "cm"),
           column_split =clini$Group)
pdf("/home/anoop/Desktop/COVID_Omics/TargettedMetabolm/HeatMap/TCA_AA.pdf")
draw(H1,merge_legend = TRUE)
dev.off()



###### FACS Correlation New gating

library(psych)
library(reshape2)

cd=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/HC_cds.txt",check.names=FALSE,header = TRUE)
facs=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/HC_facs.txt",check.names=FALSE,header = TRUE)
Corr=corr.test(as.matrix(cd),as.matrix(facs), use = "pairwise",method="spearman",adjust="none")
write.table(Corr$r,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/HC_R.txt",sep="\t",col.names = NA,quote = FALSE)
write.table(Corr$p,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/HC_P.txt",sep="\t",col.names = NA,quote = FALSE)

corR=melt(Corr$r)
corP=melt(Corr$p)
corR$Pvalue=corP$value

write.table(corR,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/HC_Results.txt",sep="\t",col.names = NA,quote = FALSE)

############

cd=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Mild_cds.txt",check.names=FALSE,header = TRUE)
facs=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Mild_facs.txt",check.names=FALSE,header = TRUE)
Corr=corr.test(as.matrix(cd),as.matrix(facs), use = "pairwise",method="spearman",adjust="none")
write.table(Corr$r,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Mild_R.txt",sep="\t",col.names = NA,quote = FALSE)
write.table(Corr$p,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Mild_P.txt",sep="\t",col.names = NA,quote = FALSE)

corR=melt(Corr$r)
corP=melt(Corr$p)
corR$Pvalue=corP$value

write.table(corR,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Mild_Results.txt",sep="\t",col.names = NA,quote = FALSE)

##########

cd=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Severe_cds.txt",check.names=FALSE,header = TRUE)
facs=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Severe_facs.txt",check.names=FALSE,header = TRUE)
Corr=corr.test(as.matrix(cd),as.matrix(facs), use = "pairwise",method="spearman",adjust="none")
write.table(Corr$r,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Severe_R.txt",sep="\t",col.names = NA,quote = FALSE)
write.table(Corr$p,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Severe_P.txt",sep="\t",col.names = NA,quote = FALSE)

corR=melt(Corr$r)
corP=melt(Corr$p)
corR$Pvalue=corP$value

write.table(corR,file="/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Severe_Results.txt",sep="\t",col.names = NA,quote = FALSE)


########


library(ComplexHeatmap)
library(circlize)
info=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Info.txt",row.names = 1)
head(info)
ha = HeatmapAnnotation(df = info,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Type=list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8),
                                                                Cohort=list(direction = "horizontal",grid_width = unit(0.5, "cm"))),
                                                      grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                       col = list(Type=c("CD4"="#ccae94","CD8"="#86959f"),Cohort=c("Mild"="#ffd700","Severe"="#ffa500","HC"="#48690E")))


col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))

pva=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Pvalues.txt",row.names = 1)
pva=as.matrix(pva)
corr=read.delim("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Correlation.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           cell_fun = function(j, i, x, y, w, h, fill) {
             gb = textGrob("*")
             gb_w = convertWidth(grobWidth(gb), "mm")
             gb_h = convertHeight(grobHeight(gb), "mm")
             if(pva[i, j] < 0.05) {
               grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
             }
           },
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),column_split = info$Cohort,row_names_side = "left",
           name = "Correlation Coeficient",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(8, "cm"),width  = unit(5, "cm"))



pdf("/home/anoop/Desktop/COVID_Omics/FeatureSelection/NEW/New_Correlation/New_Gating/Correlation.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


##################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,Condition=="HD" & CCR6 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CCR6 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Conv=ggplot()+geom_point(data=subset(data,Condition=="Conv" & CCR6 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+labs(title = "Conv")+
  geom_point(data=subset(data,Condition=="Severe" & CCR6 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Moderate=ggplot()+geom_point(data=subset(data,Condition=="Moderate" & CCR6 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CCR6 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Severe=ggplot()+geom_point(data=subset(data,Condition=="Severe" & CCR6 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+labs(title = "Severe")+
  geom_point(data=subset(data,Condition=="Severe" & CCR6 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR6),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/UMAP_CCR6.pdf",height = 13,width =15)
ggarrange(HD,Conv,Moderate,Severe,nrow = 2,ncol = 2)
dev.off()


##############

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/per.txt",header = TRUE)
head(data)

pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/bar_CCR6.pdf")
ggplot(data=data, aes(x=Cell, y=(Count),fill=Cell))+ylab("Log(Exp)")+labs(title = "CCR6")+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("CD14+ Mono"="#949494","Activated T"="#2475aa","Naive T"="#d03330",
                              "B"="#a25526","NK"="#e4e183","CD16+ Mono"="#5cb99c","GDT"="#8a4f95",
                              "MAIT"="#53a341","Platelet"="#dc86b3","Pro T"="#e87b1f","Plasma"="#bbd2e0",
                              "Mono DCs"="#4c4cff","pDC"="#919ec3","HSC"="#9fca5f"))+
  theme(axis.title = element_blank(),legend.title = element_blank(),axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=10),plot.margin = margin(1,0.5,10,0.5, "cm"),
        plot.title = element_text(hjust = 0.5,size=10),legend.position = "none")
dev.off()

#################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,Condition=="HD" & CCR5 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Conv=ggplot()+geom_point(data=subset(data,Condition=="Conv" & CCR5 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+labs(title = "Conv")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Moderate=ggplot()+geom_point(data=subset(data,Condition=="Moderate" & CCR5 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Severe=ggplot()+geom_point(data=subset(data,Condition=="Severe" & CCR5 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+labs(title = "Severe")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP_CCR5.pdf",height = 13,width =15)
ggarrange(HD,Conv,Moderate,Severe,nrow = 2,ncol = 2)
dev.off()


############

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,Condition=="HD" & CCR2 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CCR2 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Conv=ggplot()+geom_point(data=subset(data,Condition=="Conv" & CCR2 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+labs(title = "Conv")+
  geom_point(data=subset(data,Condition=="Severe" & CCR2 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Moderate=ggplot()+geom_point(data=subset(data,Condition=="Moderate" & CCR2 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CCR2 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Severe=ggplot()+geom_point(data=subset(data,Condition=="Severe" & CCR2 == 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+labs(title = "Severe")+
  geom_point(data=subset(data,Condition=="Severe" & CCR2 > 0), aes(x=UMAP_1, y=UAMP_2,color=CCR2),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP_CCR2.pdf",height = 13,width =15)
ggarrange(HD,Conv,Moderate,Severe,nrow = 2,ncol = 2)
dev.off()
###################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,Condition=="HD" & CX3CR1 == 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CX3CR1 > 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Conv=ggplot()+geom_point(data=subset(data,Condition=="Conv" & CX3CR1 == 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+labs(title = "Conv")+
  geom_point(data=subset(data,Condition=="Severe" & CX3CR1 > 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Moderate=ggplot()+geom_point(data=subset(data,Condition=="Moderate" & CX3CR1 == 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CX3CR1 > 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Severe=ggplot()+geom_point(data=subset(data,Condition=="Severe" & CX3CR1 == 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+labs(title = "Severe")+
  geom_point(data=subset(data,Condition=="Severe" & CX3CR1 > 0), aes(x=UMAP_1, y=UAMP_2,color=CX3CR1),size=0.25,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP_CX3CR1.pdf",height = 13,width =15)
ggarrange(HD,Conv,Moderate,Severe,nrow = 2,ncol = 2)
dev.off()


##############################################################

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,Condition=="HD" & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Conv=ggplot()+geom_point(data=subset(data,Condition=="Conv" & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "Conv")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Moderate=ggplot()+geom_point(data=subset(data,Condition=="Moderate" & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#f2f2f2","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Severe=ggplot()+geom_point(data=subset(data,Condition=="Severe" & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "Severe")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP_CCR5_2.pdf",height = 13,width =15)
ggarrange(HD,Conv,Moderate,Severe,nrow = 2,ncol = 2)
dev.off()


##########

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP.txt",header = TRUE)

HD=ggplot()+geom_point(data=subset(data,(Condition=="HD" & Condition=="Conv") & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))


Moderate=ggplot()+geom_point(data=subset(data,(Condition=="Moderate" & Condition=="Severe") & CCR5 < 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1,shape=19)+labs(title = "Moderate")+
  geom_point(data=subset(data,Condition=="Severe" & CCR5 >= 1), aes(x=UMAP_1, y=UAMP_2,color=CCR5),size=1.5,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#f2f2f2","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))


library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/UMAP_CCR5_3.pdf",height = 13,width =15)
ggarrange(HD,Moderate,nrow = 2,ncol = 2)
dev.off()


##########

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/CCR2/BarPlot.txt")
head(data)
M=melt(data)

pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/CCR2/BarPlot_CCR2.pdf")
ggplot(M, aes(fill=variable, y=value, x=Cell)) + ylab("% of cells expresses in cluster")+facet_wrap(~ variable)+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c(HD="#48690E",Conv="#8BB443",Moderate="#ffd700",Severe="#ffa500"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.margin = margin(2,0.5,10,0.5, "cm"),
        legend.title = element_blank(),axis.title.x = element_blank())
dev.off()  


data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/BarAll.txt")
head(data)
M=melt(data)
head(M)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/CCR6/CCr5/BarAll.pdf")
ggplot(M, aes(fill=variable, y=value, x=Cell)) + ylab("% of cells expresses in cluster")+facet_wrap(~ Gene, ncol = 1,scales ="free" )+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c(HD="#48690E",Conv="#8BB443",Moderate="#ffd700",Severe="#ffa500"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_blank(),axis.title.x = element_blank())
dev.off()  

remotes::install_version("Seurat", version = "4.0.0")

install.packages("testthat",dependencies=TRUE)

install.packages('ggplot2')

library(ggplot2)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/UMAP.txt",header = TRUE)


pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/UMAP_Severe.pdf")
ggplot(data=subset(data,Condition=="Severe"), aes(x=UMAP_1, y=UAMP_2,color=FinalCellType)) + geom_point(size=0.2,shape=19)+theme_bw()+
  scale_color_manual(values=c("CD14+ Mono"="#6497b1","Activated T"="#cc9f13","Naive T"="#fde05c",
                              "B"="#b7b81c","NK"="#b4a7d6","CD16+ Mono"="#005b96","GDT"="#466913",
                              "MAIT"="#d09c9c","Platelet"="#a6a6a6","Pro T"="#e87b1f","Plasma"="#e5f2c4",
                              "Mono DCs"="#03396c","pDC"="#800000","HSC"="#9fca5f"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=8), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),plot.margin = margin(5,2,5,4, "cm"),
                                  legend.title=element_text(size=8),legend.text=element_text(size=8),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))
dev.off()

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/UMAP_3.txt",row.names = 1)
HD=ggplot()+geom_point(data=subset(data,Condition2=="HD" & SLC25A11 == 0), aes(x=UMAP_1, y=UAMP_2),color="#f2f2f2",size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition2=="HD" & SLC25A11 > 2), aes(x=UMAP_1, y=UAMP_2,color=SLC25A11),size=0.75,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Cov=ggplot()+geom_point(data=subset(data,Condition2=="Covid" & SLC25A11 == 0), aes(x=UMAP_1, y=UAMP_2),color="#f2f2f2",size=0.25,shape=19)+labs(title = "Covid")+
  geom_point(data=subset(data,Condition2=="Covid" & SLC25A11 > 2), aes(x=UMAP_1, y=UAMP_2,color=SLC25A11),size=0.75,shape=19)+theme_bw()+
  scale_colour_gradientn(colours = c("#f2f2f2","#ff1919","#e50000","#990000","#7f0000"), breaks = c(0,1,2,3,4))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/UMAP_3.pdf",height = 13,width =15)
ggarrange(HD,Cov,nrow = 2,ncol = 1)+theme(plot.margin = margin(5,10,5,10, "cm"))
dev.off()


########

data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/UMAP_3.txt",row.names = 1)
HD=ggplot()+geom_point(data=subset(data,Condition2=="HD" & SLC25A10 == 0), aes(x=UMAP_1, y=UAMP_2),color="#cccccc",size=0.25,shape=19)+labs(title = "HD")+
  geom_point(data=subset(data,Condition2=="HD" & SLC25A10 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A10"),size=0.75,shape=19)+
  geom_point(data=subset(data,Condition2=="HD" & SLC25A1 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A1"),size=0.75,shape=19)+
  geom_point(data=subset(data,Condition2=="HD" & SLC25A11 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A11"),size=0.75,shape=19)+
  scale_color_manual(values=c(SLC25A10="#3d5404",SLC25A1="#0b3ee0",SLC25A11="#521296"))+
  theme_bw()+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

Cov=ggplot()+geom_point(data=subset(data,Condition2=="Covid" & SLC25A10 == 0), aes(x=UMAP_1, y=UAMP_2),color="#cccccc",size=0.25,shape=19)+labs(title = "Covid")+
  geom_point(data=subset(data,Condition2=="Covid" & SLC25A10 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A10"),size=0.75,shape=19)+
  geom_point(data=subset(data,Condition2=="Covid" & SLC25A1 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A1"),size=0.75,shape=19)+
  geom_point(data=subset(data,Condition2=="Covid" & SLC25A11 > 2), aes(x=UMAP_1, y=UAMP_2,color="SLC25A11"),size=0.75,shape=19)+
  scale_color_manual(values=c(SLC25A10="#3d5404",SLC25A1="#0b3ee0",SLC25A11="#521296"))+
  theme_bw()+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=15), panel.grid.major = element_blank(),axis.text = element_text(size=15),
                                  panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=15),
                                  legend.title=element_text(size=15),legend.text=element_text(size=15),legend.key.size = unit(0.5, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/UMAP_3.pdf",height = 13,width =15)
ggarrange(HD,Cov,nrow = 2,ncol = 1)+theme(plot.margin = margin(5,10,5,10, "cm"))
dev.off()


#######
library(reshape2)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/BarAll.txt")
head(data)
M=melt(data)
head(M)
pdf("/home/anoop/Desktop/COVID_Omics/Public_ScRNA/ScRNA/NEW_gene_UMAP/BarAll.pdf")
ggplot(M, aes(fill=variable, y=value, x=Cell)) + ylab("% of cells expresses in cluster")+facet_wrap(~ Gene, ncol = 1,scales ="free" )+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c(HD="#48690E",Covid="#ffa500"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title = element_blank(),axis.title.x = element_blank())
dev.off()  


###############

library(umap)
library(ggplot2)
data=read.delim("/home/anoop/Desktop/Avinash/FINAL_Analysis/Neutrophil/heatmap/figure/Normalized_data.txt",row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/COVID_Omics/AssocitationNetwork/Co-Exp/Analysis/results/DEG/Mild_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)


########## Cell type Covariates adjusted differential expression as per revision

library(DESeq2)
library(RUVSeq)

count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Severe","Healthy"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

######

count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Mild/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Mild/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Mild","Conv"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Mild/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

###
rm(list = ls())
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Severe/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Severe/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Severe","Conv"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Severe/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

##

rm(list = ls())
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Conv/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Conv/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Conv","Healthy"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Conv/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

##

rm(list = ls())
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Mild","Healthy"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)


##

rm(list = ls())
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Severe","Mild"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/DesignNew.txtResult_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

##

rm(list = ls())
count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Covid","Healthy"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)
###

count=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/DesignNew.txt",row.names = 1, check.names = FALSE)

library(dplyr)
filtered=count%>%filter(rowSums(across(where(is.numeric)))!=0)
Myset <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(Design))
Myset <- betweenLaneNormalization(Myset, which="upper")
design <- model.matrix(~Group+BMI+classical_monocyte+gdT_cell+intermediate_monocyte+MAIT_Tcell+memory_Bcell+memory_CD8_Tcell+neutrophil, data=pData(Myset))
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]

y <- DGEList(counts=counts(Myset), group=as.factor(Design$Group))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
Myset2 <- RUVr(Myset, genes, k=1, res)
dds <- DESeqDataSetFromMatrix(countData = counts(Myset2),colData = pData(Myset2),design = ~ W_1 + Group)
dds=DESeq(dds)
res=results(dds,contrast = c("Group","Covid","Healthy"),independentFiltering = FALSE)                         
write.table(as.data.frame(res),file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Result_RUV.txt",sep="\t", quote=FALSE,col.names = NA)

############## PIANO

library(piano)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Mild/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Mild/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Severe/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Conv-Severe/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Conv/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Conv/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/PianoKEGG.txt",sep="\t",col.names = NA)

##
rm(list = ls())
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/Input_Piano.txt",header = TRUE,row.names = 1)
p <- data[6]
head(p)
lfc<-data[3]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/PianoKEGG.txt",sep="\t",col.names = NA)


dat=read.delim("/home/anoop/Desktop/Slack_BackUp/test.txt",header = TRUE)
dat$BH =p.adjust(dat$p,method = "BH")
write.table(dat,file="/home/anoop/Desktop/Slack_BackUp/testBH.txt",sep="\t",col.names = NA)

?p.adjust


#############


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/WithGeneSymbol_Result_RUV.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/Volcano_HC_vs_COVID.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=2.5,label.size = NA,nudge_y = 0.5,nudge_x = 0.5,max.overlaps = 400)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1.5),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1.5, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1.5, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-7, 7.1), breaks = seq(-7, 7.1, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=10),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),plot.margin = margin(3.7,3.5,3.7,3.5, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (adj.p)")
dev.off()


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/WithGeneSymbol_Result_RUV.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/Volcano_Mild_vs_Severe.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=2.5,label.size = NA,nudge_y = 0.5,nudge_x = 0.5,max.overlaps = 600)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange <= -1.5),aes(x=log2FoldChange,y=-log10(padj )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange >= 1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange > 0 & log2FoldChange < 1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj <.05 & log2FoldChange < 0 & log2FoldChange > -1.5),aes(x=log2FoldChange,y=-log10(padj)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, padj >=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1.5, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1.5, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=10),axis.text.y=element_text(size=10),axis.text.x=element_text(size=10),plot.margin = margin(3.7,3.5,3.7,3.5, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (adj.p)")
dev.off()


################### Upset Plot

library(UpSetR)
HC_Covid=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/HC_Covid.txt",row.names = 1,header = TRUE)
HC_Conv=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/HC_Conv.txt",row.names = 1,header = TRUE)
HC_Mild=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/HC_Mild.txt",row.names = 1,header = TRUE)
HC_Severe=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/HC_Severe.txt",row.names = 1,header = TRUE)
Conv_Mild=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/Conv_Milde.txt",row.names = 1,header = TRUE)
Conv_Severe=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/Conv_Severe.txt",row.names = 1,header = TRUE)
Mild_Severe=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/Severe_Mild.txt",row.names = 1,header = TRUE)

list=list(Covid19_vs_HC=as.vector(HC_Covid$GeneName),Conv_vs_HC=as.vector(HC_Conv$GeneName),Mild_vs_HC=as.vector(HC_Mild$GeneName),Severe_vs_HC=as.vector(HC_Severe$GeneName),Mild_vs_Conv=as.vector(Conv_Mild$GeneName),Severe_vs_Conv=as.vector(Conv_Severe$GeneName),Severe_vs_Mild=as.vector(Mild_Severe$GeneName))
pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/UpsetPlot/Upset.pdf",width=10,height = 5)
upset(fromList(list),order.by = "freq",set_size.show = TRUE,set_size.scale_max=1250,main.bar.color = "#009999",sets.bar.color = "#009999",nsets = 7)
dev.off()

library(piano)
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/PianoIP.txt",header = TRUE,row.names = 2)
p <- data[7]
head(p)
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/PianoKEGG.txt",sep="\t",col.names = NA)

##############
data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/PianoHeatmap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.1,0.5,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(12, "cm"),width  = unit(3, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC_Covid/PianoHeatmap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

########


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/PianoHeatMap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.1,0.5,1,1.5,2,2.05,2.2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(12, "cm"),width  = unit(3, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Mild-Severe/PianoHeatmap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()


##########

data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MedianCellFr.txt")
M=reshape2::melt(data)
head(M)
#scale_color_manual(values=c(HC="#48690E",Mild="#ffd700",Severe="#ffa500"))+
library(ggplot2)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Permanova/MedianCellFr.txt.pdf")
ggplot(M) +
  geom_point(data=M,aes(x=variable,y=factor(ID,levels = unique(ID)),size=value,color=value))+theme_bw()+
  scale_color_gradient(low="#669cc0",high="#00365a")+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x =element_text(size=10,color="black",angle =45, vjust =1.2, hjust=1),
        axis.text.y = element_text(size=10,color="black"),plot.title = element_text(hjust = 0.5),plot.margin = margin(4,4,5,3.5, "cm"),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),
        legend.direction = "vertical",
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Median Cell Proportion",ncol = 1),color=guide_legend(title = NULL))+
  scale_size_continuous(range = c(2,5))
dev.off()


###########


Data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/BoxInput.txt",header = TRUE,check.names = FALSE)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list(c("HC", "Mild"),c("HC", "Severe"),c("Conv", "Mild"),c("Conv", "Severe"),c("Mild", "Severe"))
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/Box.pdf")
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~ variable, ncol = 3,scales ="free")+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#335b76",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons,size=0)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position =c(0.8,0.25),panel.border = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 10),
        plot.margin = margin(3.5,2,3.5,2, "cm"),strip.text.x = element_text(size = 10),
        axis.text.y = element_text(color="black",size=8),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 4),fill = guide_legend(nrow = 4))
dev.off()

Data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/BoxInput.txt",header = TRUE,check.names = FALSE)
head(Data)
my_comparisons = list(c("HC", "Mild"),c("HC", "Severe"),c("Conv", "Mild"),c("Conv", "Severe"))
P1=ggplot(Data,aes(x=factor(Group,levels = unique(Group)),y=classical_monocyte,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(title = "classical_monocyte")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())

my_comparisons = list(c("HC", "Severe"),c("Mild", "Severe"))
P2=ggplot(Data,aes(x=factor(Group,levels = unique(Group)),y=gdTcell,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(title = "gdTcell")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())

my_comparisons = list(c("Conv", "Severe"),c("HC", "Severe"),c("Mild", "Severe"))
P3=ggplot(Data,aes(x=factor(Group,levels = unique(Group)),y=MAIT_Tcell,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(title = "MAIT_Tcell")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())

my_comparisons = list(c("HC", "Severe"))
P4=ggplot(Data,aes(x=factor(Group,levels = unique(Group)),y=memory_Bcell,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(title = "memory_Bcell")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())

my_comparisons = list(c("HC", "Mild"),c("HC", "Severe"),c("Conv", "Mild"),c("Conv", "Severe"))
P5=ggplot(Data,aes(x=factor(Group,levels = unique(Group)),y=neutrophil,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, color="black",size=1,position=position_jitter(0.05))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(title = "neutrophil")+
  theme_bw()+ stat_compare_means(comparisons = my_comparisons,size=FALSE)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position =c(1.5,0.5),legend.title = element_blank(),
        axis.ticks.x=element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5,size = 9),
        panel.border = element_blank())

library(ggpubr)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/KW/NEW_KW/Box.pdf")
ggarrange(P1,P2,P3,P4,P5,nrow = 2,ncol = 3)+theme(plot.margin = margin(3.5,2,3.5,2, "cm"))
dev.off()

#########

dat=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Cluster/Coding_LCPM_UMAP.txt",row.names = 1)
gg <- merge(dat,aggregate(cbind(mean.pc1=V1,mean.pc2=V2)~Group,dat,mean),by="Group")
pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Cluster/UMAP.pdf")
ggplot(gg, aes(V1,V2,color=Group))+geom_point(size=3,aes(fill=Group))+
  geom_point(aes(x=mean.pc1,y=mean.pc2),size=5)+
  geom_segment(aes(x=mean.pc1, y=mean.pc2, xend=V1, yend=V2))+
  scale_color_manual(values=c(HC="#39540b",Conv="#6f9035",Mild="#ccac00",Severe="#cc8400"))+
  scale_fill_manual(values=c(HC="#48690E",Conv="#8BB443",Mild="#ffd700",Severe="#ffa500"))+
  labs(x="UMAP1",y="UMAP1")+
  theme(axis.title = element_text(size=10),legend.position = "bottom",plot.margin = margin(3,3,3,3, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))+
  stat_ellipse(geom = "polygon",aes(fill = Group), alpha = 0.15,linetype = "dotted")
dev.off()


##########


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/PwyHeatmap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,1,1.25,1.5,2,2.25), c("#d8d8d8" ,"#dae6f0","#7da7ca","#4682b4","#3f75a2","#386890"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Info.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.15, "cm"),
                       annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.35, "cm"),
                                                                   grid_height = unit(0.35, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                       col = list(Comparison=c("a) HC_vs_COVID19"="#555555", "b) HC_vs_Conv"="#666666","c) HC_vs_Mild"="#767676",
                                               "d) HC_vs_Severe"="#888888","e) Conv_vs_Mild"="#999999","f) Conv_vs_Severe"="#bbbbbb","g) Mild_vs_Severe"="#dddddd"),
                                  Type=c("p adj (dist.dir.dn)"="#329932", "p adj (dist.dir.up)"="#ff3232")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,column_gap = unit(0.05, "cm"),column_names_gp=gpar(fontsize = 2),
           heatmap_legend_param =list(grid_width = unit(0.35, "cm"),grid_height = unit(0.35, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 6),height  = unit(16, "cm"),width  = unit(5, "cm"),column_split = sf$Comparison)


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/PwyHeatmap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

##############################


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/All_Pwy.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.1,0.5,1,1.5,2,2.05,2.2), c("white" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/Info.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.3, "cm"),
                                                                   grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                       col = list(Comparison=c("HC_vs_COVID"="#c5cee1","Mild_vs_Severe"="#7d8daf"),
                                  Type=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 6),height  = unit(14, "cm"),width  = unit(4, "cm"),column_split = sf$Comparison)


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/All_Pwy.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

#####


library(ComplexHeatmap)
library(circlize)

sampleinfo=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HeatMap/SampleInfo.txt",row.names = 1,check.names = FALSE)

ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",
                                                                                                                                               grid_width = unit(0.75, "cm"),grid_height = unit(0.75, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("a) HC"="#48690E","b) Convalescent"="#8BB443","c) Mild"="#ffd700","d) Severe"="#ffa500")))

data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HeatMap/GeneTPM.txt",row.names = 1,check.names = FALSE)

col_fun1 = colorRamp2(c(1.5,1,0.5, 0,-0.5,-1,-1.5), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
col_fun1 = colorRamp2(c(1.5,1, 0,-1,-1.5), c("#450c5b","#453883","#21888c","#a9db35","#fbe81d"))
col_fun1 = colorRamp2(c(1.5,1,0.5, 0,-0.5,-1,-1.5), c("#071c56","#243296","#3793bb","#4cb7c5","#7bcfbd","#edf5b9","#fdfdda"))
H1=Heatmap(as.matrix(t(scale(t(data)))),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation = ha,heatmap_legend_param =list(grid_width = unit(0.75, "cm"),grid_height = unit(0.75, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(30, "cm"),width  = unit(15, "cm"),
           )

cmp=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HeatMap/comparison.txt",row.names = 1,check.names = FALSE)
ha3 = HeatmapAnnotation(df = cmp,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Comparison = list(direction = "horizontal",grid_width = unit(0.75, "cm"),grid_height = unit(0.75, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                        col = list(Comparison=c("a) HC_vs_Covid"="#dae6f0","b) HC_vs_Conv"="#b5cde1","c) HC_vs_Mild"="#7da7ca",
                                                "d) HC_vs_Severe"="#588ebb","e) Conv_vs_Mild"="#4682b4","f) Conv_vs_Severe"="#386890","g) Mild_vs_Severe"="#2a4e6c")))

col_fun_lfc = colorRamp2(c(3, 1, 0, -1, -3), c("#98120c","#fe1f14" ,"white","#9999ff","#0000ff"))
LFC=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HeatMap/LFC7.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(3, "cm"),show_row_names = FALSE,show_column_names = FALSE,
           top_annotation = ha3,heatmap_legend_param =list(grid_width = unit(0.75, "cm"),grid_height = unit(0.75, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           row_names_gp =gpar(fontsize = 5),height  = unit(45, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "white") 
t=H1+H2


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HeatMap/GeneTPM_2.pdf",height = 20,width =17)
draw(t,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


########################## Marker Gene expression heatmap
library(reshape2)
R=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/HC_R.txt",header = TRUE,sep=",")
head(R)
Cas=reshape::cast(R)
head(Cas)
write.table(Cas,file="/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/HC_casted.txt",sep="\t",col.names = NA,quote = FALSE,na = "0")


R=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/Conv_R.txt",header = TRUE)
head(R)
Cas=reshape::cast(R)
head(Cas)
write.table(Cas,file="/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/Conv_casted.txt",sep="\t",col.names = NA,quote = FALSE,na = "0")

##################
library(ComplexHeatmap)
library(circlize)
com=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/HC_Type.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Type = list(direction = "horizontal",grid_width = unit(0.3, "cm"),
                                                              grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 0), 
                                                              labels_gp = gpar(fontsize = 8))),
                   col = list(Type=c("CD4+_memory_T_cell"="#26c9f8","Classical_monocyte"="#0197fc",
                                     "Eosinophil"="#e4c100","gdT-cell"="#9a9cfa","Intermediate_monocyte"="#2d3af6",
                                     "MAIT_T-cells"="#ffcb65","Memory_B_cell"="#f90567","Memory_CD8+_T-cell"="#fa9798",
                                     "Myeloid_dendritic_cell"="#f19f81","Naive_B_cell"="#ce6503","Naive_CD4+_T_cell"="#f76639",
                                     "Naive_CD8+_T_cell"="#B9e8BA","Natural_killer_cell"="#9f019d","Neutrophil"="#cdcb00",
                                     "Non-classical_monocyte"="#0303d9","Plasmacytoid_dendritic_cell"="#359590","Regulatory_T_(Treg)_cell"="#ffc0e9")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Type = list(direction = "horizontal",grid_width = unit(0.3, "cm"),
                                                                   grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 0), 
                                                                   labels_gp = gpar(fontsize = 8))),
                        col = list(Type=c("CD4+_memory_T_cell"="#26c9f8","Classical_monocyte"="#0197fc",
                                          "Eosinophil"="#e4c100","gdT-cell"="#9a9cfa","Intermediate_monocyte"="#2d3af6",
                                          "MAIT_T-cells"="#ffcb65","Memory_B_cell"="#f90567","Memory_CD8+_T-cell"="#fa9798",
                                          "Myeloid_dendritic_cell"="#f19f81","Naive_B_cell"="#ce6503","Naive_CD4+_T_cell"="#f76639",
                                          "Naive_CD8+_T_cell"="#B9e8BA","Natural_killer_cell"="#9f019d","Neutrophil"="#cdcb00",
                                          "Non-classical_monocyte"="#0303d9","Plasmacytoid_dendritic_cell"="#359590","Regulatory_T_(Treg)_cell"="#ffc0e9")))




corr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/HC_casted.txt",row.names = 1)
head(corr)
col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#e50000","#ff0000" ,"white","#6666ff", "#0000ff"))
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),column_dend_height = unit(0.5, "cm"),
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha3,left_annotation = ha,
           row_split = com$Cell,column_split = com$Cell,
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(8, "cm"))

head(corr)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/HC.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

############


com=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/Conv_Type.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Type = list(direction = "horizontal",grid_width = unit(0.3, "cm"),
                                                              grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 0), 
                                                              labels_gp = gpar(fontsize = 8))),
                   col = list(Type=c("CD4+_memory_T_cell"="#26c9f8","Classical_monocyte"="#0197fc",
                                     "Eosinophil"="#e4c100","gdT-cell"="#9a9cfa","Intermediate_monocyte"="#2d3af6",
                                     "MAIT_T-cells"="#ffcb65","Memory_B_cell"="#f90567","Memory_CD8+_T-cell"="#fa9798",
                                     "Myeloid_dendritic_cell"="#f19f81","Naive_B_cell"="#ce6503","Naive_CD4+_T_cell"="#f76639",
                                     "Naive_CD8+_T_cell"="#B9e8BA","Natural_killer_cell"="#9f019d","Neutrophil"="#cdcb00",
                                     "Non-classical_monocyte"="#0303d9","Plasmacytoid_dendritic_cell"="#359590","Regulatory_T_(Treg)_cell"="#ffc0e9")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Type = list(direction = "horizontal",grid_width = unit(0.3, "cm"),
                                                                   grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 0), 
                                                                   labels_gp = gpar(fontsize = 8))),
                        col = list(Type=c("CD4+_memory_T_cell"="#26c9f8","Classical_monocyte"="#0197fc",
                                          "Eosinophil"="#e4c100","gdT-cell"="#9a9cfa","Intermediate_monocyte"="#2d3af6",
                                          "MAIT_T-cells"="#ffcb65","Memory_B_cell"="#f90567","Memory_CD8+_T-cell"="#fa9798",
                                          "Myeloid_dendritic_cell"="#f19f81","Naive_B_cell"="#ce6503","Naive_CD4+_T_cell"="#f76639",
                                          "Naive_CD8+_T_cell"="#B9e8BA","Natural_killer_cell"="#9f019d","Neutrophil"="#cdcb00",
                                          "Non-classical_monocyte"="#0303d9","Plasmacytoid_dendritic_cell"="#359590","Regulatory_T_(Treg)_cell"="#ffc0e9")))




corr=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/Conv_casted.txt",row.names = 1)
head(corr)
col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#e50000","#ff0000" ,"white","#6666ff", "#0000ff"))
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),column_dend_height = unit(0.5, "cm"),
           row_dend_width = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha3,left_annotation = ha,
           row_split = com$Cell,column_split = com$Cell,
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(8, "cm"))

head(corr)
pdf("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/NewHeatMaps/Convalescent.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


######

data=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/PWy/PWY_Count.txt")
head(data)
MM=melt(data)
head(MM)

pdf("/home/anoop/Desktop/COVID_Omics/EPIC/Network/HC_and_Conv/PWy/PWY_Count.pdf",width = 12,height = 8)
ggplot(MM, aes(fill=variable, y=value, x=Term)) + geom_bar(position='dodge', stat='identity')+ylab("Number of Genes present")+
  scale_fill_manual(values=c(HC="#48690E",Convalescent="#8BB443",Mild="#FFD700",Severe="#FFA500"))+theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1,size=10),
                                                                                  plot.margin = margin(3,0.5,6,2, "cm"),axis.title.x = element_blank())
dev.off()



#########


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/pwy_HeatMap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.1,0.5,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(8, "cm"),width  = unit(2.5, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Mild/pwy_HeatMap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

####


data=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/PwyHeatmap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.1,0.5,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))

sf=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/HC_Conv-Mild/Label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (dist.dir.dn)"="#004C00",
                                          "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (mix.dir.up)"="#E50000","e) p adj (dist.dir.up)"="#7F0000")))

H1=Heatmap(as.matrix(-log10(data)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),top_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 7),height  = unit(12, "cm"),width  = unit(2.5, "cm"))


pdf("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq_Adjusted/HC-Severe/PwyHeatmap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()
