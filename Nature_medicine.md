Single-cell immune landscape of human atherosclerotic plaques
https://doi.org/10.1038/s41591-019-0590-4

```R
list.files("./data_used/")
[1] "merged_plaque_gex-mnn_macrophages.txt.gz"              "merged_plaque_gex-mnn_t_cells.txt.gz"                 
[3] "merged_plaque_gex-umi-data-mnn-cat_macrophages.txt.gz" "merged_plaque_gex-umi-data-mnn-cat_t_cells.txt.gz" 

#文件批量读取，然后构建Seurat对象后保存为Rdata 文件，供后续使用
rm(list=ls())
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(clustree)) 
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))
library(harmony)
# 1-----
path = "data_used/" 
fs=list.files(path = 'data_used/',pattern="gex-mnn")

sceList = lapply(fs, function(f){
  gex <- read.table(file.path('data_used',f),sep="\t",header=T,row.names=1)
  list <- c()
  for (i in colnames(gex)){
       a <- strsplit(i,"[.]")[[1]][2]
       list <- append(list,a)
  }
  colnames(gex) <- list
    
  sce = CreateSeuratObject(counts = gex)
  sce$group=strsplit(f,'_')[[1]][2]
  sce$celltype=gsub('.txt.gz','',strsplit(f,'_')[[1]][4])
  sce
})
#合并
sce_all <- merge(sceList[[1]], y= sceList[-1] )
#备份数据
save(sce_all,file="sce_all.Rdata")

# 2.预处理数据 ------

macro <- read.table("./data_used/merged_plaque_gex-mnn_macrophages.txt.gz",sep="\t",header=T,row.names=1)
    celltype1 <- c()
    sample1 <- c()
    Symptomatic1 <- c()
for (i in colnames(macro)){
   a <- paste(strsplit(i,'[..]')[[1]][16],strsplit(i,'[..]')[[1]][17],sep="")
   celltype1 <- append(celltype1,a)
   m <- strsplit(i,"[.]")[[1]][7] 
   sample1 <- append(sample1,m)
   n <- strsplit(i,"[.]")[[111]
   Symptomatic1 <- append(Symptomatic1,n)
}

t <- read.table("./data_used/merged_plaque_gex-mnn_t_cells.txt.gz",sep="\t",header=T,row.names=1)
      celltype2 <- c()
      sample2 <- c()
      Symptomatic2 <- c()
for (i in colnames(t)){
   a <- paste(strsplit(i,'[..]')[[1]][16],strsplit(i,'[..]')[[1]][17],sep="")
   celltype2 <- append(celltype2,a)
   m <- strsplit(i,"[.]")[[1]][7] 
   sample2 <- append(sample2,m)
   n <- strsplit(i,"[.]")[[1]][11]
   Symptomatic2 <- append(Symptomatic2,n)
}
celltype <- c(celltype1,celltype2)
sample <- c(sample1,sample2)
Symptomatic <- c(Symptomatic1,Symptomatic2)
sce_all$subcelltype <- celltype
sce_all$sample <- sample
sce_all$Symptomatic <- Symptomatic


sce=sce_all
sce <- NormalizeData(sce)
sce = FindVariableFeatures(sce)
sce = ScaleData(sce, features = rownames(sce))
# 3.降维分析-----
# 3.1 PCA降维 ------
sce = RunPCA(sce, npcs = 50)


sce <- RunHarmony(sce,group.by.vars="sample",plot_convergence=TRUE)
sce <- RunUMAP(sce,reduction="harmony",dims=1:20)
sce <- RunTSNE(sce,reduction="harmony",dims=1:20)

DimPlot(sce,reduction="umap",label.size=5,label=T,repel=T,pt.size=1)
DimPlot(sce,reduction="tsne",label.size=5,label=T,repel=T,pt.size=1)
                          
                          
                          
                
                          
                          
                          
                          

#3.2 umap和tsne降维 ------
sce = RunTSNE(sce, dims = 1:20)
sce = RunUMAP(sce, dims = 1:20)
# 4. 细胞聚类分析 -----
# 4.1 聚类
sce=FindNeighbors(sce, dims = 1:20, prune.SNN = 1/15)
#clustering 分析（设置不同的分辨率）
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce=FindClusters(sce, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
# 4.2 可视化 ------
# 4.2.1 树 ------
(p1_tree=clustree(sce@meta.data, prefix = "RNA_snn_res."))

# 4.2.2 res = 1 umap可视化 ------
sce <- SetIdent(sce, value = "RNA_snn_res.1") 
#UMAP可视化
(p2 <- DimPlot(sce , reduction = "umap",label = T, label.box = T) + NoAxes())
(p3 <- DimPlot(sce, reduction = "umap", group.by = "celltype", 
               label = T, label.box = T) +  NoAxes())

#5 注释-----
genes_to_check <- rownames(sce)
ct = as.data.frame(sce@assays$RNA@counts)
(p4 <- DotPlot(sce, col.min=0, dot.min = 0.9,
              features = unique(genes_to_check),assay='RNA') + coord_flip())

#命名
celltype=data.frame(ClusterID=0:12,celltype='unkown')
celltype[celltype$ClusterID == 3,2]='macrophages'
celltype[celltype$ClusterID %in% c(2,4,5,6,8,9,10,11),2]='CD4+'
celltype[celltype$ClusterID %in% c(0,7,12),2]='CD8+' 
celltype[celltype$ClusterID == 1,2]='CD4+CD8+' 
#将celltype信息添加值meta.data
sce@meta.data$CELLTYPE = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'CELLTYPE'] <- celltype$celltype[i]
}
# 细胞注释结果umap可视化 ------
(p5 <- DimPlot(sce, reduction = "umap", group.by = "CELLTYPE", 
               label = T, label.box = T) +  NoAxes())
```
