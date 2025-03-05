getwd()

setwd("/home/OSUMC.EDU/li176/work/Jianying/Jocye")

library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
library(uwot)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(data.table)
library(enrichR)
library(tidyverse)
library(limma)
library(harmony)
library(cowplot)
library(here)
library(qs)
library(Polychrome)

dat.one<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/RA/HRR059379")
dat.two<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/RA/HRR059382")
dat.three<-Read10X(data.dir ="/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/RA/HRR059385")
dat.four<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/RA/HRR059388")
dat.five<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/Healthy Control/HRR059375")
dat.six<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/Healthy Control/HRR059376")
dat.seven<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/Healthy Control/HRR059377")
dat.eight<-Read10X(data.dir = "/fs/ess/PAS2177/Jianying/Fianal/Chinese/Data/Healthy Control/HRR059378")


objone<-CreateSeuratObject(counts = dat.one,project = "RA",min.cells = 3,min.features = 200)
objtwo<-CreateSeuratObject(counts = dat.two,project = "RA",min.cells = 3,min.features = 200)
objthree<-CreateSeuratObject(counts = dat.three,project = "RA",min.cells = 3,min.features = 200)
objfour<-CreateSeuratObject(counts = dat.four,project="RA",min.cells = 3, min.features = 200)
objfive<-CreateSeuratObject(counts = dat.five,project = "RA",min.cells = 3, min.features = 200)
objsix<-CreateSeuratObject(counts = dat.six,project = "RA",min.cells = 3,min.features = 200)
objseven<-CreateSeuratObject(counts = dat.seven,project = "RA",min.cells = 3,min.features = 200)
objeight<-CreateSeuratObject(counts = dat.eight,project = "RA",min.cells = 3,min.features = 200)

objone$Condition="RA"
objtwo$Condition="RA"
objthree$Condition="RA"
objfour$Condition="RA"
objfive$Condition="HC"
objsix$Condition="HC"
objseven$Condition="HC"
objeight$Condition="HC"

objone$Tissue="PBMC"
objtwo$Tissue="PBMC"
objthree$Tissue="PBMC"
objfour$Tissue="PBMC"
objfive$Tissue="PBMC"
objsix$Tissue="PBMC"
objseven$Tissue="PBMC"
objeight$Tissue="PBMC"

objone$Serotype="positive"
objtwo$Serotype="positive"
objthree$Serotype="positive"
objfour$Serotype="positive"
objfive$Serotype="positive"
objsix$Serotype="positive"
objseven$Serotype="positive"
objeight$Serotype="positive"



objone[["percent.mt"]]<-PercentageFeatureSet(objone,pattern = "^MT-")
objone<-subset(objone,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objone<-NormalizeData(objone, normalization.method = "LogNormalize", scale.factor = 10000)
objone<-FindVariableFeatures(objone,selection.method = "vst",nfeatures = 2000)

objtwo[["percent.mt"]]<-PercentageFeatureSet(objtwo,pattern = "^MT-")
objtwo<-subset(objtwo,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objtwo<-NormalizeData(objtwo, normalization.method = "LogNormalize", scale.factor = 10000)
objtwo<-FindVariableFeatures(objtwo,selection.method = "vst",nfeatures = 2000)

objthree[["percent.mt"]]<-PercentageFeatureSet(objthree,pattern = "^MT-")
objthree<-subset(objthree,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objthree<-NormalizeData(objthree, normalization.method = "LogNormalize", scale.factor = 10000)
objthree<-FindVariableFeatures(objthree,selection.method = "vst",nfeatures = 2000)

objfour[["percent.mt"]]<-PercentageFeatureSet(objfour,pattern = "^MT-")
objfour<-subset(objfour,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objfour<-NormalizeData(objfour, normalization.method = "LogNormalize", scale.factor = 10000)
objfour<-FindVariableFeatures(objfour,selection.method = "vst",nfeatures = 2000)

objfive[["percent.mt"]]<-PercentageFeatureSet(objfive,pattern = "^MT-")
objfive<-subset(objfive,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objfive<-NormalizeData(objfive, normalization.method = "LogNormalize", scale.factor = 10000)
objfive<-FindVariableFeatures(objfive,selection.method = "vst",nfeatures = 2000)

objsix[["percent.mt"]]<-PercentageFeatureSet(objsix,pattern = "^MT-")
objsix<-subset(objsix,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objsix<-NormalizeData(objsix, normalization.method = "LogNormalize", scale.factor = 10000)
objsix<-FindVariableFeatures(objsix,selection.method = "vst",nfeatures = 2000)

objseven[["percent.mt"]]<-PercentageFeatureSet(objseven,pattern = "^MT-")
objseven<-subset(objseven,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objseven<-NormalizeData(objseven, normalization.method = "LogNormalize", scale.factor = 10000)
objseven<-FindVariableFeatures(objseven,selection.method = "vst",nfeatures = 2000)

objeight[["percent.mt"]]<-PercentageFeatureSet(objeight,pattern = "^MT-")
objeight<-subset(objeight,subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
objeight<-NormalizeData(objeight, normalization.method = "LogNormalize", scale.factor = 10000)
objeight<-FindVariableFeatures(objeight,selection.method = "vst",nfeatures = 2000)


obj.list<-list(objone, objtwo, objthree,objfour,objfive,objsix,objseven,objeight)
Tdat <-
    RunHarmony(Tdat, "obj.list", plot_convergence = FALSE)
Tdat <-
    RunUMAP(Tdat, reduction = "harmony", dims = 1:20)
Tdat <-
    FindNeighbors(Tdat, reduction = "harmony", dims = 1:20)
Tdat <- FindClusters(Tdat, resolution = 0.2)



data.anchors<-FindIntegrationAnchors(obj.list, dims = 1:20)
M.dat<-IntegrateData(anchorset = data.anchors, dims = 1:20)
# correct the gene name
counts <- GetAssayData(M.dat, assay = "RNA",slot = "counts")
#delete the GRCH38- in the gene label use gsub function
rownames(counts) <- gsub(pattern="GRCh38-", replacement = "" ,x=rownames(counts))
rownames(counts)
#creat a new seurat object
Tdat <- CreateSeuratObject(counts = counts,meta.data = M.dat@meta.data)

table(Tdat$Tissue)
table(Tdat$Serotype)
table(Tdat$Condition)

condition.list <- SplitObject(Tdat, split.by = "Condition")
#vnormalize and identify variable features for each dataset independently
condition.list<- lapply(X = condition.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = condition.list)
immune.anchors <- FindIntegrationAnchors(object.list = condition.list, anchor.features = features,dims =1:30)
# intergrate the data
Tdat <-IntegrateData(anchorset = immune.anchors,dims =1:30)
DefaultAssay(Tdat) <- "integrated"
Tdat <- ScaleData(Tdat, verbose = FALSE)
Tdat <- RunPCA(Tdat, npcs = 30, verbose = FALSE)
Tdat <- RunUMAP(Tdat, reduction = "pca", dims = 1:20)
Tdat <- FindNeighbors(Tdat, reduction = "pca", dims = 1:30, k.param = 35)
Tdat <- FindClusters(Tdat, resolution = 0.3)


#################################use Harmony to integrate the data#######################
#############Harmony#################################################################
########################################################################################
Idents(Tdat) <- Tdat$Condition
DefaultAssay(Tdat) <- "RNA"
selected_mito <- WhichCells(Tdat, expression = percent_mito < 0.2)
Tdat <- ScaleData(Tdat, verbose = FALSE)
Tdat <-
    FindVariableFeatures(Tdat,
                         selection.method = "vst",
                         nfeatures = 2000)
Tdat <-RunPCA(
    Tdat,
    pc.genes = Tdat@var.genes,
    npcs = 20,
    verbose = FALSE
)
Tdat <-
    RunHarmony(Tdat, "Condition", plot_convergence = FALSE)
Tdat <-
    RunUMAP(Tdat, reduction = "harmony", dims = 1:20)
Tdat <-
    FindNeighbors(Tdat, reduction = "harmony", dims = 1:20)
Tdat <- FindClusters(Tdat, resolution = 0.2)

DefaultAssay(Tdat) <- "RNA"
Idents(Tdat) <- Tdat$seurat_clusters

DimPlot(
    Tdat,
    reduction = "umap",
    group.by = "seurat_clusters",
    split.by = "orig.ident",
    label = T
)

####################################################################






# relabeld cluster 
# Rename Cluster


#############################################
# create Tfh name
save(dat.Tfh,file = "dat.Tfh.RData")
save(Tdat,file = "Tdat.RData")
DefaultAssay(Tdat) <- "RNA"
pc <- WhichCells(Tdat, expression =BCL6>0 & CD4> 0) #& BCL6>0)
#pc <- WhichCells(Tdat, expression = CXCR5 >0 & PDCD1>0)
Tdat$Tfh<- ifelse(colnames(Tdat) %in% pc, "Tfh", "other")
Idents(Tdat) <- "Tfh"

Condition_Tfh<- paste(Tdat$Condition, "_",Tdat$Tfh)
names(Condition_Tfh) <- colnames(x = Tdat)
Tdat <- AddMetaData(
    object = Tdat,
    metadata = Condition_Tfh,
    col.name = 'Condition_Tfh')

Idents(Tdat)<-"Condition_Tfh"
table(Tdat$Condition_Tfh)
##############subset Tfh################

dat.Tfh <- subset(Tdat, idents=c("HC _ Tfh","RA _ Tfh"))

genenew.list<-c("ITGB7","RORC","ICOS","MAF","SH2D1A","SLAMF1","LY9","CD84","SLAMF6","IL2RA","IL2","IL21")

png(filename = "Human_Tfh_CD4+BCL6+.png", width = 9, height = 6, res = 300, units = "in")
DotPlot(dat.Tfh,features = PP4, group.by="Condition")+theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
dev.off()
############################ load the data


PP4 <- fread("Th17derTfh_signgenesPP.csv")
PP4<- PP4[,2]
PP4



SP4 <- fread("Th17derTfh_signgenesSpleen.csv")
SP4<- SP4[,2]
SP4
## USE The PP1 and SP1#################################
DefaultAssay(object = dat.Tfh) <- "RNA"
dat.Tfh<-AddModuleScore(object = dat.Tfh, features = PP4, name = "PPTfhhuman_score")
head(dat.Tfh)
viodat<-data.frame(Condition=dat.Tfh$Condition, PPTfhhuman_score= dat.Tfh$PPTfhhuman_score1)
fun_mean<-function(x){return(data.frame(y=mean(x),label=round(mean(x,na.rm=T),3)))}

compare_means(PPTfhhuman_score~Condition, data = viodat)

# my_comparisons<-list(c("HC _ CD4-", "HC _ CD4+"),c("HC _ CD4-", "RA _ CD4-"),c("HC _ CD4-", "RA _ CD4+"),c("RA _ CD4+", "HC _ CD4+" ), c("RA _ CD4-", "HC _ CD4+"),c("RA _ CD4-", "RA _ CD4+"))

#pdf(filename = "PP-2 signature gene in CD4+BCL6+ Tfh cell.pdf", width = 4, height = 4, units = "in")
table(viodat$Condition)
table(viodat$PPTfhhuman_score)

p2 <- ggplot(viodat, aes(x=Condition,y=PPTfhhuman_score, fill=Condition))+
    geom_violin()+
    scale_fill_manual(values = c("#56B4E9","#FF9999"))+
    geom_boxplot(width=0.1,fill="white")+
    stat_compare_means(comparisons = my_comparisons,label="p.signif")+
    stat_compare_means(label.y = 3.0, label.x = 1.5)+
    ylim(-1, 3)+
    stat_summary(fun.data = fun_mean,geom = "text",vjust=-0.1)+
    theme_classic()+labs(title = "PP-2 signature gene in CD4+BCL6+ Tfh cell")+theme(axis.text = element_text(size = 25))+theme(axis.title = element_text(size = 25))+
    theme(legend.text = element_text(size=25)) +theme(legend.title = element_text(size=25))
p2
#dev.off()


DefaultAssay(dat.Tfh) <- "RNA"
dat.Tfh<-AddModuleScore(object = dat.Tfh, features = SP4, name = "sTfh_score")
head(dat.Tfh)
viodat3<-data.frame(Condition=dat.Tfh$Condition, sTfh_score= dat.Tfh$sTfh_score1)
fun_mean<-function(x){return(data.frame(y=mean(x),label=round(mean(x,na.rm=T),3)))}

compare_means(sTfh_score~Condition, data = viodat1)
#my_comparisons<-list(c("HC _ CD4-", "HC _ CD4+"),c("HC _ CD4-", "RA _ CD4-"),c("HC _ CD4-", "RA _ CD4+"),c("RA _ CD4+", "HC _ CD4+" ), c("RA _ CD4-", "HC _ CD4+"),c("RA _ CD4-", "RA _ CD4+"))


#pdf(filename = "spleen-2 signature gene in CD4+BCL6+  Tfh cell.pdf", width = 6, height = 8, units = "in")
p1 <- ggplot(viodat3, aes(x=Condition,y=sTfh_score, fill=Condition))+
    geom_violin()+
    scale_fill_manual(values = c("#56B4E9","#FF9999"))+
    geom_boxplot(width=0.1,fill="white")+
    stat_compare_means(comparisons = my_comparisons,label="p.signif")+
    stat_compare_means(label.y = 3.0, label.x = 1.5)+
    ylim(-1, 3.0)+
    stat_summary(fun.data = fun_mean,geom = "text",vjust=0.5)+
    theme_classic()+
    labs(title = "spleen-2 signature gene in  CD4+BCL6+ Tfh cell")+theme(axis.text = element_text(size = 25))+theme(axis.title = element_text(size = 25))+
    theme(legend.text = element_text(size=25))+theme(legend.title = element_text(size=25))  
p1
#dev.off() 
#####################################################
#pdf(filename = "spleen and PP signature gene in CD4+BCL6+  Tfh cell.pdf", width = 6, height = 8, units = "in")
p1+p2
#dev.off()

#

