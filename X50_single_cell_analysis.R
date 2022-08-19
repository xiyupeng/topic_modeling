### 
# Author: Xiyu Peng
# single cell analysis for flow cytometry data
#
# Codes below requires at least 250G memory space 
#


library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(future)
library(tidyverse)

start_time <- Sys.time()
Sys.time()

### 4 cores 
plan("multiprocess", workers = 80)
plan()
options(future.globals.maxSize = 100000 * 1024^2)


# the working directory
setwd("/gpfs/mskmind_ess/pengx1/flow/X50_rerun")
set.seed(42) ## reproducible 

mapping_file <-"/gpfs/mskmind_ess/pengx1/flow/analysis/17-162_file_mapping.csv"
meta<-read.csv(mapping_file)
meta$rds_suffix<-paste(meta$SAMPLE_ID,"_intensity_dat_preGated.rds",sep = '')

submeta<-meta 
count_matrix<-NULL
downsample<-30000000
remove_negative<-0

### pool all samples together
for(i in 1:nrow(submeta)){
  raw_file<-submeta$rds_suffix[i]
  data_submat<-NULL
  tmp<-try(data_submat<-readRDS(raw_file),silent = TRUE)
  if(!is.null(data_submat)){
    data_submat<-as.data.frame(data_submat)
    data_submat<-data_submat[data_submat$cd3_pos == 1,7:34] ##  only keep cd3+ cells and useful columns
    if(nrow(data_submat) > downsample){
      data_submat<-data_submat[sample(rownames(data_submat),size = downsample, replace = F),]
    }
    if(remove_negative){
      data_submat[data_submat<0]<-NA
      data_submat<-na.omit(data_submat)
    }
    data_submat$time<-submeta$TP[i]
    data_submat$patient<-submeta$PATIENT[i]
    data_submat$sample<-submeta$SAMPLE_ID[i]
    count_matrix<-rbind(count_matrix,data_submat)
    rm(data_submat)
  }
}

data_mat_cd3pos<-count_matrix
rm(count_matrix)

setwd("/gpfs/mskmind_ess/pengx1/flow/analysis")

data_real_marks<-as.matrix(data_mat_cd3pos[,-c(10,29:31)])
dim(data_real_marks)
colnames(data_real_marks)


## relable T cells (needed for creating seurat object)
new_id<-as.character(1:nrow(data_real_marks))
rownames(data_real_marks)<-new_id

####---------------------------------------------------------------------------------
## Analysis

#### create seurat object 

pool_X50 <- CreateSeuratObject(counts = t(data_real_marks), project = "17-162-test")
pool_X50
pool_X50[['time']]<-factor(data_mat_cd3pos$time)
pool_X50[['patient']]<-factor(data_mat_cd3pos$patient)
pool_X50[['sample']]<-factor(data_mat_cd3pos$sample)

### quality control (remove following samples, due to lack of cells)
all.samples<-levels(pool_X50$sample)
rm_samples<-c("17-162-32C","17-162-34B","17-162-07B")
rm_samples
all.samples<-all.samples[! all.samples %in% rm_samples]
pool_X50<-subset(pool_X50,subset = sample %in% all.samples)

## remove object to save memory
rm(data_real_marks)
rm(data_mat_cd3pos)

## select genes (may think about to remove the step)
pool_X50 <- FindVariableFeatures(object = pool_X50)
VariableFeatures(pool_X50)

### scaled data
all.genes <- rownames(pool_X50)
pool_X50 <- ScaleData(pool_X50, features = all.genes)

### linear dimension reduction 
pool_X50 <- RunPCA(pool_X50, features = all.genes,approx=FALSE,npcs = 26)
DimPlot(pool_X50, reduction = "pca")
ElbowPlot(pool_X50,ndims = 26)

Sys.time()

## umap for visualization (based on knn graph)
pool_X50 <- RunUMAP(pool_X50, dims = 1:26,min.dist = 0.1)

Sys.time()


### do graph-based clustering based on pca result
pool_X50 <- FindNeighbors(pool_X50, dims = 1:26, nn.eps = 0.5,k.param = 5)

Sys.time()

### Run clustering under different resolution
res<-c(0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3)

pool_X50 <- FindClusters(pool_X50, resolution = res)

colnames(pool_X50@meta.data)

## choose resolution = 1.5
Idents(pool_X50)<-"RNA_snn_res.1.5"

### use visualization to check
#DimPlot(pool_X50,reduction = "umap", label = TRUE)
#DimPlot(pool_X50,reduction = "umap",split.by = "time")
#DimPlot(pool_X50,reduction = "umap",split.by = "patient")

### save the seurat object
saveRDS(pool_X50, file = "seurat_object_17162_rerun_all_dim26.rds")

#proc.time()
end_time <- Sys.time()
end_time - start_time

### codes below are for further analysis

### UMAP plot
DimPlot(subset(pool_X50, idents = 0:19),reduction = "umap", label = TRUE)

## Feature plots
FeaturePlot(subset(pool_X50, idents = 0:19), features = c("CD4","CD8","FOXP3","CD45RA","CCR7","LAG3","PD1","CTLA4","KI67"), ncol =3, min.cutoff = 0, slot = "scale.data",max.cutoff = 3)  
FeaturePlot(subset(pool_X50, idents = 0:19), features = c("CD57","EOMES","TIGIT","TBET","GZM-B","CXCR5","GITR","CD38","CD27"), ncol =3, min.cutoff = 0, slot = "scale.data",max.cutoff = 3)
FeaturePlot(subset(pool_X50, idents = 0:19), features = c("CD28","CD127","ICOS","CD25","CCR4","TIM3","HLADR"), ncol =3, min.cutoff = 0, slot = "scale.data",max.cutoff = 3)


### save the marker average expression matrix for heatmap
genes.reorder<-c("CD4","CD8","FOXP3","CCR7","CD45RA","GITR","CD25","CCR4","CD27","CD28","CD127","TBET","GZM-B","CD57","EOMES","LAG3","CXCR5","HLADR","KI67","TIM3","CD38","ICOS","CTLA4","PD1","TIGIT")
subset_object<-subset(pool_X50, downsample = 10000)
all.genes<-rownames(subset_object)
mat1<-NULL
for (gene in all.genes){
  mat1<-rbind(mat1,tapply(subset_object[['RNA']]@scale.data[gene,],Idents(subset_object),mean))
}
rownames(mat1)<-all.genes
save(mat1,file = "average_cluster.Rdata")


### save seurat object of given patients
patients<-c("17-162-27","17-162-05","17-162EXT-09")

for(pt in patients){

    pt_seurat<-subset(pool_X50,subset = patient == pt)   
    pt_seurat$pool_cluster<-pt_seurat$RNA_snn_res.1.5  ##keep previous clustering result
    pt_seurat <- RunUMAP(pt_seurat, dims = 1:26,min.dist = 0.3)
    saveRDS(pt_seurat,file= paste0("/gpfs/mskmind_ess/pengx1/flow/analysis/seurat_pt",pt,"redo.rds"))
}

### ridge plots

RidgePlot(object = pool_X50, features = 'KI67',idents = c(4,8,9,12,16),sort = "increasing",cols = c('4' = "#C09B00" , '8' = "#00BB4E", '9' = "#00BF7D", '12' = "#00BAE0",'16' = "#C77CFF" ))


sessionInfo()

