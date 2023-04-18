##    title: "SCRGOT 2023 Coder Upgrade Session 4: Combining single cell data sets"
##    author: "Xin Wang"
##    date: '2023-5-3'
##    email: xin.wang@nationwidechildrens.org
##    output: html_document

## Method Part1: Integrating using the canonical correlation analysis (CCA) in Seruat
## Method reference: https://satijalab.org/seurat/articles/integration_introduction.html

## Dataset Description:
## This tutorial has four samples, 
## including two replicates from human normal kidney (GSM5627690, GSM5627691), 
## two replicates from  autosomal dominant polycystic kidney disease (GSM5627695 and GSM5627696).
## The dataset can be downloaded from GEO:GSE185948 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185948
## To reduce the time, We randomly select 1000 of cells in each sample.
## Dataset reference: Muto Y, Dixon EE, Yoshimura Y, Wu H et al. Defining cellular complexity in human autosomal dominant polycystic kidney disease by multimodal single cell analysis. Nat Commun 2022 Oct 30;13(1):6497. PMID: 36310237

## please download the datasets from the github 
## github datasets: 


## Step 1: Load the packages, 
## Please install all the packages before loading them.
library(Seurat)
library(gridExtra)
library(ggplot2)
library(tidyverse)

## Step 2: Set up the workspace environment:
setwd("/Users/XXW004/Documents/Projects/LectureForSingleCellGroup/IntegrationSection/TestDatasets/")

## check the files in the current directory
Datafiles<-list.files(path = "./", recursive = F, full.names = F)
Datafiles

## Step 3: Read each sample
Control1<- readRDS(file = "GSM5627690_cont1.small.rds")
Control2<- readRDS(file = "GSM5627691_cont2.small.rds")
Disease1<- readRDS(file = "GSM5627695_PKD1.small.rds")
Disease2<- readRDS(file = "GSM5627696_PKD2.small.rds")

## briefly check wether we successfully load the single cell files
dim(Control1@meta.data)
str(Control1@meta.data)
head(Control1@meta.data)

## Step 4: Merge these objects
## merge multiple objects from the list, we also set up the cell id 
scrna<-merge(x=Control1, y=c(Control2, Disease1,Disease2), add.cell.ids = c("A","B","C","D"), project="Integration")
scrna@meta.data

## we then create a group column fro the meta data based on their condition
## mutate to add a group into the meta, 
## str_split split the cell ID and select the first element, simplify =TRUE will return a character matrix

scrna@meta.data <- scrna@meta.data %>% mutate(group = stringr::str_split(row.names(scrna@meta.data), "_", simplify = TRUE)[,1])

## please check the string split function:
test<-str_split(row.names(scrna@meta.data), "_", simplify = TRUE)[,1]
tail(scrna@meta.data)
# check the structure of meta data
str(scrna@meta.data)
levels(factor(scrna@meta.data$group))
# view the meta data
View(scrna@meta.data)

## Step 5: QC & filtering

# calculate mitochondrial percentage

scrna$mitoPercent <-PercentageFeatureSet(scrna, pattern = '^MT-')

# let's check the quality of datasets by eveluating the mitochondrial percentage, number of Counts, number of features.
head(scrna@meta.data)
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA"))
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA") , split.by = 'group')

# filtering
scrna <- subset (scrna, subset =mitoPercent <10 & nFeature_RNA >500 & nCount_RNA >200 )

## Step 6: Perform the standard workflow
# perform the standard workflow to figure out if there are any batch effects
scrna<- NormalizeData(object = scrna)
scrna<- FindVariableFeatures(object = scrna)
scrna<- ScaleData(object = scrna)
scrna<- RunPCA(object = scrna)
ElbowPlot(scrna)
scrna<- FindNeighbors(object = scrna, dim=1:15)
scrna<- FindClusters(object=scrna)
scrna<- RunUMAP(object = scrna, dims = 1:15)


## Step 6: Plot the merged UMAP
p1 <- DimPlot(scrna, reduction = 'umap', group.by = 'group')
p1

p2<- DimPlot(scrna, reduction = 'umap', split.by = 'group')

grid.arrange(p1, p2, ncol=2, nrow=1)

# Obviously, there are batch effects for the sample, so we perform integration to correct for batch effects

## Step 6: Split the objects, normalization, and find varaible feature for each objects
# split the objects
SplitedObjects<- SplitObject(scrna, split.by = 'group')

# check the split objects
SplitedObjects
length(SplitedObjects)

## Normalized dataset and Find variable features
for (i in 1: length(SplitedObjects)){
  SplitedObjects[[i]] <-NormalizeData(object = SplitedObjects[[i]])
  SplitedObjects[[i]] <- FindVariableFeatures(object = SplitedObjects[[i]])
}

## Step7: Prepare the integration feature, integration achor
# select integration features

features<- SelectIntegrationFeatures(object.list = SplitedObjects)
head(features)

# find integration anchor (CCA)

anchors<- FindIntegrationAnchors(object.list = SplitedObjects, anchor.features = features)

## Step 8: Run the integration using CCA
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
# Plot PCA
PCAPlot(seurat.integrated, split.by= "group")

# we can see from the PCA that a good overlay of several condtions by PCA

# Now, we can also visualize with UMAP.

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15, reduction = "pca")

## Step 9: Visualize the UMAP
# Plot UMAMP
DimPlot(seurat.integrated)

#Idents(seurat.integrated)<-ordered(factor(Idents(seurat.integrated)), levels=c(0:16))

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'group')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'group',
              cols = c('red','green','blue','yellow'))
#p5 <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'group',
 #             cols = c('red','green','blue','yellow'))

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow=2)

## side by side spitted by groups

DimPlot(seurat.integrated, reduction = 'umap', split.by = 'group')

