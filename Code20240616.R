library(dplyr)
library(tidyr)
library(raster)
library(dplyr)
library(gplots)
library(data.table)
library(ggplot2)
library(GenVisR)
library(descriptr)
library(vioplot)
library(pheatmap)
library(survival)
library(corrplot)
library(preprocessCore)
source("Cibersort.R")
library(e1071)
library(tidyr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggimage)
library(plyr)
library(ggord)
library(yyplot)
library("sva")
library(WGCNA)

GPL14951_anno=data.table::fread("GPL14951-11332.txt",skip = 28) #璇诲彇娉ㄩ噴锛宻kip琛ㄧず璺宠繃寮�澶寸殑
GPL14951_anno=GPL14951_anno[,c(1,12)]
GSE89632_matrix=read.table("GSE89632_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE89632_matrix$ID=rownames(GSE89632_matrix)
GSE89632=merge(x=GSE89632_matrix,y=GPL14951_anno,by="ID",all.x = T)
GSE89632$ID=NULL
rowMeans=apply(GSE89632,1,function(x)mean(as.numeric(x),na.rm=T))
GSE89632=GSE89632[order(rowMeans,decreasing = T),]
GSE89632=GSE89632[!duplicated(GSE89632[,dim(GSE89632)[2]]),]
GSE89632=GSE89632[!grepl("///",GSE89632$Symbol),]
rownames(GSE89632)=GSE89632[,dim(GSE89632)[2]]
GSE89632=GSE89632[,-dim(GSE89632)[2]]
write.csv(GSE89632,"GSE89632_Symbol_expr.csv",quote=F)

GPL11532_anno=data.table::fread("GPL11532-32230.txt",skip = 12) #璇诲彇娉ㄩ噴锛宻kip琛ㄧず璺宠繃寮�澶寸殑
GPL11532_anno=GPL11532_anno[,c(1,10)]
GPL11532_anno=separate(data = GPL11532_anno, col = gene_assignment, into = c("ensemble", "Symbol","function1"), sep = "//",convert = T)
GPL11532_anno=GPL11532_anno[,-c(2,4)]
GSE48452_matrix=read.table("GSE48452_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE48452_matrix$ID=rownames(GSE48452_matrix)
GSE48452=merge(x=GSE48452_matrix,y=GPL11532_anno,by="ID",all.x = T)
GSE48452$ID=NULL
GSE48452=na.omit(GSE48452)
GSE48452$Symbol=trim(GSE48452$Symbol)
rowMeans=apply(GSE48452,1,function(x)mean(as.numeric(x),na.rm=T))
GSE48452=GSE48452[order(rowMeans,decreasing = T),]
GSE48452=GSE48452[!duplicated(GSE48452[,dim(GSE48452)[2]]),]
GSE48452=GSE48452[!grepl("///",GSE48452$Symbol),]
rownames(GSE48452)=GSE48452[,dim(GSE48452)[2]]
GSE48452=GSE48452[,-dim(GSE48452)[2]]
write.csv(GSE48452,"GSE48452_Symbol_expr.csv",quote=F)

GSE89632$ID=rownames(GSE89632)
GSE48452$ID=rownames(GSE48452)
pca1=merge(x=GSE48452,y=GSE89632,by="ID",all.x = T)
rownames(pca1)=pca1$ID
pca1$ID=NULL
write.csv(pca1,"pca1.csv",quote=F)
pca1=na.omit(pca1)
#ppca1=cor(pca1,method = "pearson") 
pca_df=as.data.frame(t(pca1))
pca_result1=prcomp(pca_df,center=T,scale=F)
exprfl=read.csv("DEGfenlei.csv",header=T,sep = ",",row.names = 1,stringsAsFactors = F)
ggord(pca_result1,grp_in=exprfl$type,arrow=0,vec_ext=0,txt=NULL,cols=c("red","blue"),size=1,
      ellipse=T,poly = F)



pca2=as.matrix(pca1)
csif=read.csv("csif.csv",header = T,sep=",",row.names = 1)
modcombat = model.matrix(~1, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=pca2, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata, "hb_jz.csv", sep = ",", quote = F)

pca_df2=as.data.frame(t(combat_edata))
pca_result2=prcomp(pca_df2,center=T,scale=F)
exprfl=read.csv("DEGfenlei.csv",header=T,sep = ",",row.names = 1,stringsAsFactors = F)
ggord(pca_result2,grp_in=exprfl$type,arrow=0,vec_ext=0,txt=NULL,cols=c("red","blue"),size=1,
      ellipse=T,poly = F)

data=as.data.frame(combat_edata)
write.table(data, "data.csv", sep = ",", quote = F)

#宸紓鍒嗘瀽
library(limma)
fenzu=read.csv("group.csv")
data=data[,fenzu$sample]
datahcss=data[,c(1:72)]

group1=read.csv("group1.csv",sep = ",",row.names = 1,header = T,check.names = F)
design=model.matrix(~0+factor(group1$type)) #type璁剧疆鎴愪竴涓猰odel 鐭╅樀
colnames(design)=levels(factor(group1$type))
rownames(design)=colnames(datahcss)
fit=lmFit(datahcss,design) #绾挎�ф嫙鍚?
cont.matrix=makeContrasts(Steatosis-Control,levels = design)
fit2=contrasts.fit(fit,cont.matrix)  #鐢ㄥ姣旀ā鍨嬭繘琛屽樊鍊艰绠?
fit2=eBayes(fit2) #璐濆彾鏂楠?
temp0output=topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by = "B",resort.by = "M")
temp0output$ID=rownames(temp0output)

hdb=read.csv("autophagy.csv",header = T)
HDB_diffhcss=temp0output[c(temp0output$ID%in%hdb$gene),]

write.csv(HDB_diffhcss,"HDB_diffhcss.csv")

write.csv(temp0output,"limmaOut.csv")
foldchange=0.5
pvalue=0.05
temp0output1=temp0output
diffhcss=temp0output1
diffhcss=diffhcss[(diffhcss$P.Value<pvalue&(diffhcss$logFC>foldchange | diffhcss$logFC<(-foldchange))),]
write.csv(datahcss,"datahcss.csv")


x=temp0output
x$label=rownames(x)
logFCcut=0.5
pvaluecut=0.05
x[,7]=ifelse((x$P.Value<pvaluecut & x$logFC>logFCcut),"red",ifelse((x$P.Value<0.05 & x$logFC<(-logFCcut)),"blue","grey30"))
size=ifelse((x$P.Value<pvaluecut & abs(x$logFC)>logFCcut),4,2)
xmin=-1
xmax=1
ymax=15
ymin=0
valo=ggplot(data=x,aes(x=logFC,y=-log10(P.Value),label=label)) +
  geom_point(alpha=0.6,size=size,colour=x[,7]) +
  scale_color_manual(values=c("lightgrey","navy","red")) +
  labs(x=bquote(~log[2]~"(fold Change)"),y=bquote(~-log[10]~italic("p-value"))) +
  ylim(c(ymin,ymax)) +
  scale_x_continuous(
    breaks = c(-1.5,-1,-logFCcut,0,logFCcut,1,1.5),
    labels = c(-1.5,-1,-logFCcut,0,logFCcut,1,1.5),
    limits = c(-2,2)
  ) +
  geom_vline(xintercept = logFCcut,color="grey40",linetype="longdash",size=0.5) +
  geom_vline(xintercept = -logFCcut,color="grey40",linetype="longdash",size=0.5) +
  geom_hline(yintercept = -log10(pvaluecut),color="grey40",linetype="longdash",size=0.5) +
  guides(colour=guide_legend(override.aes = list(shape=16))) +
  theme_bw(base_size=12,base_family = "Times") +
  theme(legend.position = "right",
        panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(face="bold",color="black",family = "Times_New_Roman",size=8),
        plot.title=element_text(hjust=0.8),
        axis.text.x=element_text(face="bold",color="black",size=15),
        axis.text.y=element_text(face="bold",color="black",size=15))

#鐑浘 
library(pheatmap)
p1=datahcss[rownames(diffhcss),]
pheatmap(p1,cluster_cols = F,cluster_rows=F,scale="row",color = colorRampPalette(c("navy","white","firebrick3"))(100),
         main = "DEGheatmap",annotation_col = group1,show_colnames = F,show_rownames = F)

library(WGCNA)
GSE89632$ID=NULL
WGCNA=GSE89632
WGCNA$CV=apply(WGCNA,1,function(x) sd(x)/mean(x))
WGCNA=WGCNA[(WGCNA$CV>=0.05),]
WGCNA$CV=NULL
WGCNA=as.data.frame(t(WGCNA))
sampleTree = hclust(dist(WGCNA), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
WGCNA=WGCNA[rownames(WGCNA)!=c("GSM2385767"),]
WGCNA=WGCNA[-62,]
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(WGCNA, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#鏋勫缓缃戠粶
cor <- WGCNA::cor
net = blockwiseModules(WGCNA, power =15 ,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3,)
cor<-stats::cor 
table(net$colors)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
q=as.data.frame(mergedColors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];


#鐩稿叧鍥?
nGenes = ncol(WGCNA);
nSamples = nrow(WGCNA);
MEs0 = moduleEigengenes(WGCNA, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,group_WGCNA, use = "p");
group_WGCNA=read.csv("group_wgcna.csv",sep=",",header = T,row.names = 1)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(9, 10, 3, 3));
# Display the correlation values within a heatmap 
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(group_WGCNA),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

geneModuleMembership = as.data.frame(cor(WGCNA, MEs, use = "p"));

brown_gene=names(WGCNA)[moduleColors=="brown"]
brown_gene=as.data.frame(brown_gene)
write.csv(brown_gene,"brown_gene.csv")

blue_gene=names(WGCNA)[moduleColors=="blue"]
blue_gene=as.data.frame(blue_gene)
write.csv(blue_gene,"blue_gene.csv")

nafld.activity.score=as.data.frame(group_WGCNA$nafld.activity.score)

geneTraitSignificance_nc=as.data.frame(cor(WGCNA, nafld.activity.score, use = "p"))

modNames = substring(names(MEs),3)
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;


#鐢籋ub鍩哄洜鏁ｇ偣鍥?
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance_nc[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

hub_UC<- abs(geneModuleMembership$MEblue)>0.8 & abs(geneTraitSignificance_nc)>0.2
hub_UC=as.data.frame(hub_UC)
hub_UC$blue_gene=rownames(hub_UC)
hub_UC=merge(hub_UC,blue_gene,by="blue_gene",all.x=F)
rownames(hub_UC)=hub_UC$blue_gene
hub_UC=hub_UC[hub_UC$`group_WGCNA$nafld.activity.score`==TRUE,]
write.csv(hub_UC,"hub_UC.csv")

write.table(GSE89632,"GSE89632.txt",quote = F,sep="\t")
cs_result=CIBERSORT("LM22.txt","GSE89632.txt",perm = 100,QN=T)
cs_result=as.data.frame(cs_result)
cs_result=cs_result[fenzu$sample[c(15:38,53:72)],]
write.csv(cs_result,"cs_result.csv")
cs_result=read.csv("cs_result.csv",header = T,row.names = 1)
cs_result=cs_result[(cs_result$P.value<0.05),] #绛涢�塒鍊?<0.05
cs_result_cl=cs_result[,-c(23,24,25)] #鍘绘帀鍚庝笁鍒?
cs_result_t=t(cs_result_cl)


#绠辩嚎鍥炬瘮杈?
library(tidyverse)
cs_result_xx=cs_result_cl%>%rownames_to_column("sample")
cs_result_xx$Group=NULL
cs_result_xx$Group=c(rep("HC",24),rep("SS",20))
library(ggpubr)
library(tidyr)
library(ggsci)
library(ggplot2)
cs_result_plot=gather(cs_result_xx,key = CIBERSORT,value=Proportion,-c(Group,sample))
ggboxplot(cs_result_plot,x="CIBERSORT",y="Proportion",
          fill="Group",palette = "Lancet") +
  stat_compare_means(aes(group=Group),
                     method = "t.test",
                     label="p.signif",
                     symnum.args = list(cutpoints=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*"," "))) +
  theme(text=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust = 1))


#鍩哄洜鍜屽厤鐤粏鑳?
sig_gene <- c("BAG1","CXCR4","PPP1R15A","CCL2","MYC","FOS")
library(psych)
xi=as.data.frame(t(GSE89632))
xi=xi[fenzu$sample[c(15:38,53:72)],]
xi <- xi[,sig_gene]
yi <- cs_result_cl

di <- corr.test(xi,yi,use="complete",method = 'spearman')

ri <- di$r
pi <- di$p

library(ggcorrplot)
ggcorrplot(t(di$r), show.legend = T, 
           p.mat = t(di$p.adj), digits = 1,  sig.level = 1,insig = 'blank',lab = T)


GPL570_anno=data.table::fread("GPL570-55999.txt",skip = 16)
GPL570_anno=GPL570_anno[,c(1,11)] 
GSE109597_matrix=read.table("GSE109597_series_matrix.txt",header = T,sep="\t",row.names = 1,comment.char = "!")
GSE109597_matrix$ID=rownames(GSE109597_matrix)
GSE109597=merge(x=GSE109597_matrix,y=GPL570_anno,by="ID",all.x = T)
GSE109597$ID=NULL
rowMeans=apply(GSE109597,1,function(x)mean(as.numeric(x),na.rm=T))
GSE109597=GSE109597[order(rowMeans,decreasing = T),]
GSE109597=GSE109597[!duplicated(GSE109597[,dim(GSE109597)[2]]),]
GSE109597=GSE109597[!grepl("///",GSE109597$`Gene Symbol`),]
rownames(GSE109597)=GSE109597[,dim(GSE109597)[2]]
GSE109597=GSE109597[,-dim(GSE109597)[2]]


group2=read.csv("group2.csv",sep = ",",row.names = 1,header = T,check.names = F)

GSE109597=GSE109597[,c(rownames(group2))]

write.csv(GSE109597,"GSE109597.csv")
design=model.matrix(~0+factor(group2$type)) #type璁剧疆鎴愪竴涓猰odel 鐭╅樀
colnames(design)=levels(factor(group2$type))
rownames(design)=colnames(GSE109597)
fit=lmFit(GSE109597,design) #绾挎�ф嫙鍚?
cont.matrix=makeContrasts(OB-Control,levels = design)
fit2=contrasts.fit(fit,cont.matrix)  #鐢ㄥ姣旀ā鍨嬭繘琛屽樊鍊艰绠?
fit2=eBayes(fit2) #璐濆彾鏂楠?
temp0output=topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by = "B",resort.by = "M")

temp0output$ID=rownames(temp0output)

HDB_GSE109597=temp0output[c(temp0output$ID%in%HDB$gene),]

write.csv(HDB_GSE109597,"HDB_GSE109597.csv")


HDB_liver=read.csv("HDB_diffhcss.csv",header = T,row.names = 1)

HDB_liver=HDB_liver[c(rownames(HDB_GSE109597)),]

write.csv(HDB_liver,"HDB_liver123.csv")


GSE213621=read.table("GSE213621TPM.txt",header = T,row.names=1)
clinc=read.csv("clincdataGSE213621.csv",header = T)
GSE213621=GSE213621[,colnames(GSE213621)%in%clinc$Sample[1:299]]
Consens_gene=GSE213621[c("BAG1","CXCR4","FOS","MYC","CCL2","PPP1R15A"),]
Consens_gene = sweep(Consens_gene,1, apply(Consens_gene,1,median,na.rm=T)) #中值中心化
Consens_gene=as.matrix(Consens_gene)

set.seed(123)
Consens_result <- ConsensusClusterPlus(Consens_gene, maxK = 10, reps = 1000, clusterAlg = "hc", distance = "pearson")
Julei=Consens_result[[2]][['consensusClass']]
icl <- calcICL(Consens_result)
icl[["clusterConsensus"]] 

write.csv(Julei,"JuleiGSE213621.csv")
write.csv(GSE213621,"GSE213621genesymbol.csv")

icl <- calcICL(Consens_result)
icl[["clusterConsensus"]]

tpm_data <- read.csv("GSE213621genesymbol.csv", row.names = 1)

group_info <- read.csv("JuleiGSE213621.csv", header = TRUE, stringsAsFactors = FALSE)


boxplot(tpm_data,outline=FALSE, notch=T, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(tpm_data)
boxplot(exprSet,outline=FALSE, notch=T,las=2)
#判断数据是否需要转换

ex=exprSet
qx=as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

group=read.csv("JuleiGSE213621.csv",sep = ",",row.names = 1,header = T,check.names = F)
exprSet=as.data.frame(exprSet)
exprSet=exprSet[,rownames(group)]

design=model.matrix(~0+factor(group$Group)) #type设置成一个model 矩阵
colnames(design)=levels(factor(group$Group))
rownames(design)=colnames(exprSet)
fit=lmFit(exprSet,design) #线性拟合
cont.matrix=makeContrasts(type2-type1,levels = design)
fit2=contrasts.fit(fit,cont.matrix)  #用对比模型进行差值计算
fit2=eBayes(fit2) #贝叶斯检验
limmaout=topTable(fit2,coef=1,n=Inf,adjust="BH",sort.by = "B",resort.by = "M")

write.csv(design,"design.csv")

foldchange=0.5
pvalue=0.05
diff=limmaout
diff=diff[(diff$P.Value<pvalue&(diff$logFC>foldchange | diff$logFC<(-foldchange))),]

write.csv(diff,"diff.csv")


task=tpm_data[rownames(diff),]
task=task[,rownames(group)]
task=as.data.frame(t(task))
write.csv(task,"task.csv")


library(cmapR)
gct_object <- parse_gctx("H:/NAFLD/camP/my_analysis.sig_queryl1k_tool.652ddb87330a010014eff389/ncs.gct")
gct_object@rdesc
camp=as.data.frame(gct_object@rdesc)
gct_object@mat


HC1=read.csv("GSM3714747_chow1_filtered_gene_bc_matrices.csv")
HC2=read.csv("GSM3714748_chow2_filtered_gene_bc_matrices.csv")
HC=merge(HC1,HC2,by="X")
write.csv(HC,"HC.csv")

a=read.csv("HC.csv",header=T,row.names = 1)
rownames(a)=a$X
a$X=NULL
library(Seurat)
library(celldex)
library(scRNAseq)
library(SingleR)
library(reshape)
library(reshape2)
pbmc=CreateSeuratObject(counts = a,min.cells = 100,min.features = 50)
a=NULL

pbmc[["precent.mt"]]=PercentageFeatureSet(pbmc,pattern="^MT-")
VlnPlot(pbmc,features = c("nFeature_RNA","nCount_RNA","precent.mt"))
pbmc=NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
pbmc=FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 3000)
all.genes=rownames(pbmc)
pbmc=ScaleData(pbmc,features = all.genes)
pbmc=RunPCA(pbmc,features = VariableFeatures(object=pbmc))
print(pbmc[["pca"]], dims = 1:30, nfeatures = 5)
pbmc=FindNeighbors(pbmc,dims=1:20)
pbmc=FindClusters(pbmc,resolution = 0.2)
pbmc=RunTSNE(object=pbmc,dims = 1:20)
DimPlot(pbmc,reduction = "tsne")

hpca.se=celldex::BlueprintEncodeData()

table(hpca.se@colData@listData$label.fine)
cluster_pbmc=pbmc@meta.data$seurat_clusters
table(pbmc@meta.data$seurat_clusters)
datas=as.SingleCellExperiment(pbmc)
pred.hesc=SingleR(test=datas,ref=hpca.se,
                  labels=hpca.se$label.main,clusters=cluster_pbmc,
                  assay.type.test="logcounts",assay.type.ref="logcounts")
table(pred.hesc$labels)
celltype=data.frame(cluster=rownames(pred.hesc),celltype=pred.hesc$labels)

for(i in 1:12177){
  index=pbmc@meta.data$seurat_clusters[i]
  pbmc@meta.data$celltype[i]=celltype[index,2]
}




pbmc.markers=FindAllMarkers(pbmc,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
write.table(pbmc.markers,file="pbmc.marker.txt",sep="\t",row.names = F,quote = F)
library(tidyverse)
sig.markers=pbmc.markers%>%select(gene,everything())%>%
  subset(p_val<0.05&abs(pbmc.markers$avg_log2FC)>1)
for(i in 1:nrow(sig.markers)){
  clusterID=as.numeric(sig.markers[i,]$cluster)
  celltypes=celltype[clusterID,2]
  sig.markers$celltype[i]=celltypes
}
write.csv(sig.markers,"sig.markers.csv")


library(viridis)
VlnPlot(object=pbmc,features = c("Myc","Ccl2","Cxcr4","Fos","Ppp1r15a","Bag1"),
        group.by = "seurat_clusters")

VlnPlot(object=pbmc,features = c("Myc","Ccl2","Cxcr4","Fos","Ppp1r15a","Bag1"),
        group.by = "celltype")

FeaturePlot(object = pbmc,features = c("Myc","Ccl2","Cxcr4","Fos","Ppp1r15a","Bag1"),cols = cividis(10))

saveRDS(pbmc,"pbmc.rds")

library(AUCell)
cells_rankings <- AUCell_buildRankings(pbmc@assays$RNA@data)
h <- read.csv("自噬gmt.csv")
geneSets<-lapply(unique(h$term),function(x){h$Genes[h$term==x]})
names(geneSets) <- unique(h$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 1, aucMaxRank=nrow(cells_rankings)*0.1)
length(rownames(cells_AUC@assays@data$AUC))
geneSet <- "Autophagy_pathway"  
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])  #提取这个通路在每一个细胞的得分
pbmc$AUC <- aucs
df<- data.frame(pbmc@meta.data, pbmc@reductions$tsne@cell.embeddings)

class_avg <- df %>%
  group_by(celltype) %>%        #这里可以改成cluster  seurat_clusters/或者其他的annotation
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
ggplot(df, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="H") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()

?AUCell_exploreThresholds
cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist=TRUE, assign=TRUE,thrP = 0.05)
geneSetName <- rownames(cells_AUC)[grep("Pyroptosis_pathway", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.05)
abline(v=0.05)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>0.05))
length(newSelectedCells)


df1=df[rownames(df)%in%newSelectedCells,]

class_avg <- df1 %>%
  group_by(celltype) %>%        #这里可以改成cluster  seurat_clusters/或者其他的annotation
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
ggplot(df1, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="H") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()


df2=df
df2$color=ifelse((df2$AUC>0.6),"red",ifelse((df2$AUC<0.6),"grey30","grey30"))

size=ifelse((df$AUC>0.6),4,2)
class_avg <- df2 %>%
  group_by(celltype) %>%        #这里可以改成cluster  seurat_clusters/或者其他的annotation
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
ggplot(df2, aes(tSNE_1, tSNE_2))  +
  geom_point(alpha=0.6,size=size,colour=df2[,11]) +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()

df4=df[rownames(df)%in%cells_assignment$Pyroptosis_pathway$assignment,]
table(pbmc@meta.data$celltype)





save(cells_rankings, file="cells_rankings.RData")

