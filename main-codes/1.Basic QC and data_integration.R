
### load required packages 
library(Seurat)
library(dplyr)
library(cowplot)
library(DoubletFinder)
library(ggplot2)

### read inputdata 
E45.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_species/scRNA_pig/cellranger_out/E45/filtered_feature_bc_matrix")
E55.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_species/scRNA_pig/cellranger_out/E55/filtered_feature_bc_matrix")
E65.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_species/scRNA_pig/cellranger_out/E65/filtered_feature_bc_matrix")
E75.data <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_species/scRNA_pig/cellranger_out/E75/filtered_feature_bc_matrix")

colnames(x = E45.data) <- paste('E45', colnames(x = E45.data), sep = '_')
colnames(x = E55.data) <- paste('E55', colnames(x = E55.data), sep = '_')
colnames(x = E65.data) <- paste('E65', colnames(x = E65.data), sep = '_')
colnames(x = E75.data) <- paste('E75', colnames(x = E75.data), sep = '_')

## Construct Seurat object
E45 <- CreateSeuratObject(counts = E45.data, project = "E45", min.cells = 3, min.features = 200)
E45[["percent.mt"]]<-PercentageFeatureSet(E45,pattern = "^MT")
VlnPlot(E45, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E45 <- subset(E45, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 0.5)
E45 <- NormalizeData(E45, normalization.method = "LogNormalize", scale.factor = 10000)
E45 <- FindVariableFeatures(E45, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E45$group<-"E45"
all.genes <- rownames(E45)
E45 <- ScaleData(E45, features = all.genes)
E45@meta.data$time <- "E45"


E55 <- CreateSeuratObject(counts = E55.data, project = "E55", min.cells = 3, min.features = 200)
E55[["percent.mt"]]<-PercentageFeatureSet(E55,pattern = "^MT")
VlnPlot(E55, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E55 <- subset(E55, subset = nFeature_RNA > 1000 &  nFeature_RNA < 6000 & percent.mt < 0.5)
E55 <- NormalizeData(E55, normalization.method = "LogNormalize", scale.factor = 10000)
E55 <- FindVariableFeatures(E55, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E55$group<-"E55"
all.genes <- rownames(E55)
E55 <- ScaleData(E55, features = all.genes)
E55@meta.data$time <- "E55"


E65 <- CreateSeuratObject(counts = E65.data, project = "E65", min.cells = 3, min.features = 200)
E65[["percent.mt"]]<-PercentageFeatureSet(E65,pattern = "^MT")
VlnPlot(E65, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E65 <- subset(E65, subset = nFeature_RNA > 1000 &  nFeature_RNA < 5000 & percent.mt < 0.5)
E65 <- NormalizeData(E65, normalization.method = "LogNormalize", scale.factor = 10000)
E65 <- FindVariableFeatures(E65, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E65$group<-"E65"
all.genes <- rownames(E65)
E65 <- ScaleData(E65, features = all.genes)
E65@meta.data$time <- "E65"

E75 <- CreateSeuratObject(counts = E75.data, project = "E75", min.cells = 3, min.features = 200)
E75[["percent.mt"]]<-PercentageFeatureSet(E75,pattern = "^MT")
VlnPlot(E75, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0) ## 5x3
E75 <- subset(E75, subset = nFeature_RNA > 1000 &  nFeature_RNA < 6000 & percent.mt < 0.5)
E75 <- NormalizeData(E75, normalization.method = "LogNormalize", scale.factor = 10000)
E75 <- FindVariableFeatures(E75, selection.method = "vst", nfeatures = 2000) ## 2000 features in default
E75$group<-"E75"
all.genes <- rownames(E75)
E75 <- ScaleData(E75, features = all.genes)
E75@meta.data$time <- "E75"

#################################### remove douplets cells with Doublefinder #################################
################################ https://www.jianshu.com/p/b1947c4156ad   ####################################
#########################   https://github.com/chris-mcginnis-ucsf/DoubletFinder #############################

############################################### Module 1 ##################################################
E45<-RunPCA(E45)
E45<-RunUMAP(E45,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_E45 <- paramSweep_v3(E45, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E45)
sweep.stats_E45 <- summarizeSweep(sweep.res.list_E45, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E45)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_E45 <- paramSweep_v3(E45, PCs = 1:20, sct = FALSE)
# gt.calls <- E45@meta.data[rownames(sweep.res.list_E45[[1]]), "GT"]
# sweep.stats_E45 <- summarizeSweep(sweep.res.list_E45, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON11 <- find.pK(sweep.stats_E45)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
E45<-FindNeighbors(E45,reduction ="pca",dims = 1:20)
E45<-FindClusters(E45,resolution = 0.5)
head(E45@meta.data)
annotationscon_E45<-E45@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E45)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(E45@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E45 <- doubletFinder_v3(E45, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_con_E45@meta.data)
# seu_con_E11@meta.data$DF.classifications_0.25_0.09_1413
table(seu_con_E45$DF.classifications_0.25_0.09_356)
# Doublet Singlet 
#  356    4397 

seu_con_E45@meta.data$cellfilter <- seu_con_E45@meta.data$DF.classifications_0.25_0.09_356
seu_con_E45@meta.data <-seu_con_E45@meta.data[,-10]

# Doublet Singlet 
# 523    6453
seu_con_E45@meta.data$time<- "E45"

DimPlot(seu_con_E45, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E45")
        
############################################### Module 1 ##################################################


############################################### Module 2 ##################################################
E55<-RunPCA(E55)
E55<-RunUMAP(E55,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_E55 <- paramSweep_v3(E55, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E55)
sweep.stats_E55 <- summarizeSweep(sweep.res.list_E55, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E55)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_E55 <- paramSweep_v3(E55, PCs = 1:20, sct = FALSE)
# gt.calls <- E55@meta.data[rownames(sweep.res.list_E55[[1]]), "GT"]
# sweep.stats_E55 <- summarizeSweep(sweep.res.list_E55, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON11 <- find.pK(sweep.stats_E55)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
E55<-FindNeighbors(E55,reduction ="pca",dims = 1:20)
E55<-FindClusters(E55,resolution = 0.5)
head(E55@meta.data)
annotationscon_E55<-E55@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E55)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(E55@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E55 <- doubletFinder_v3(E55, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_con_E55@meta.data)
# seu_con_E11@meta.data$DF.classifications_0.25_0.09_1413
table(seu_con_E55$DF.classifications_0.25_0.09_399)
# Doublet Singlet 
# 399    4919

seu_con_E55@meta.data$cellfilter <- seu_con_E55@meta.data$DF.classifications_0.25_0.09_399
seu_con_E55@meta.data <-seu_con_E55@meta.data[,-10]

# Doublet Singlet 
# 523    6453
seu_con_E55@meta.data$time<- "E55"

DimPlot(seu_con_E55, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E55")
############################################### Module 2 ##################################################


############################################### Module 3 ##################################################
E65<-RunPCA(E65)
E65<-RunUMAP(E65,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_E65 <- paramSweep_v3(E65, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E65)
sweep.stats_E65 <- summarizeSweep(sweep.res.list_E65, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E65)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_E65 <- paramSweep_v3(E65, PCs = 1:20, sct = FALSE)
# gt.calls <- E65@meta.data[rownames(sweep.res.list_E65[[1]]), "GT"]
# sweep.stats_E65 <- summarizeSweep(sweep.res.list_E65, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON11 <- find.pK(sweep.stats_E65)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
E65<-FindNeighbors(E65,reduction ="pca",dims = 1:20)
E65<-FindClusters(E65,resolution = 0.5)
head(E65@meta.data)
annotationscon_E65<-E65@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E65)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(E65@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E65 <- doubletFinder_v3(E65, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_con_E65@meta.data)
# seu_con_E11@meta.data$DF.classifications_0.25_0.09_1413
table(seu_con_E65$DF.classifications_0.25_0.09_282)
# Doublet Singlet 
# 282    3480

seu_con_E65@meta.data$cellfilter <- seu_con_E65@meta.data$DF.classifications_0.25_0.09_282
seu_con_E65@meta.data <-seu_con_E65@meta.data[,-10]

# Doublet Singlet 
# 523    6453
seu_con_E65@meta.data$time<- "E65"

DimPlot(seu_con_E65, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E65")
############################################### Module 3 ##################################################

############################################### Module 4 ##################################################
E75<-RunPCA(E75)
E75<-RunUMAP(E75,dims = 1:10)
## pK Identification (no ground-truth)-----------------------------------------------------------------------------------------
sweep.res.list_E75 <- paramSweep_v3(E75, PCs = 1:20, sct = FALSE)
head(sweep.res.list_E75)
sweep.stats_E75 <- summarizeSweep(sweep.res.list_E75, GT = FALSE)
bcmvn_CON11 <- find.pK(sweep.stats_E75)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_E75 <- paramSweep_v3(E75, PCs = 1:20, sct = FALSE)
# gt.calls <- E75@meta.data[rownames(sweep.res.list_E75[[1]]), "GT"]
# sweep.stats_E75 <- summarizeSweep(sweep.res.list_E75, GT = TRUE, GT.calls = gt.calls)
# bcmvn_CON11 <- find.pK(sweep.stats_E75)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
E75<-FindNeighbors(E75,reduction ="pca",dims = 1:20)
E75<-FindClusters(E75,resolution = 0.5)
head(E75@meta.data)
annotationscon_E75<-E75@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotationscon_E75)   ## ex: annotations <- SeuratOb@meta.data$ClusteringResults
nExp_poi <- round(0.075*(length(E75@active.ident)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_con_E75 <- doubletFinder_v3(E75, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(seu_con_E75@meta.data)
# seu_con_E11@meta.data$DF.classifications_0.25_0.09_1413
table(seu_con_E75$DF.classifications_0.25_0.09_500)
#Doublet Singlet 
#500    6160

seu_con_E75@meta.data$cellfilter <- seu_con_E75@meta.data$DF.classifications_0.25_0.09_500
seu_con_E75@meta.data <-seu_con_E75@meta.data[,-10]

seu_con_E75@meta.data$time<- "E75"

DimPlot(seu_con_E75, group.by = "cellfilter",
        cols = c("#92D713","#832688"))  + ggtitle("E75")
############################################### Module 4 ##################################################



###################################### I need to remove the douplets  #######################################
#########################################  Extract the singlets  ############################################
head(seu_con_E45@meta.data)
E45.singlet <- subset(seu_con_E45, subset = cellfilter == 'Singlet')
table(E45.singlet@meta.data$cellfilter)   ##  4397

head(seu_con_E55@meta.data)
E55.singlet <- subset(seu_con_E55, subset = cellfilter == 'Singlet')
table(E55.singlet@meta.data$cellfilter)   ##  4947

head(seu_con_E65@meta.data)
E65.singlet <- subset(seu_con_E65, subset = cellfilter == 'Singlet')
table(E65.singlet@meta.data$cellfilter)   ##  3480 

head(seu_con_E75@meta.data)
E75.singlet <- subset(seu_con_E75, subset = cellfilter == 'Singlet')
table(E75.singlet@meta.data$cellfilter)   ##  6160 

## Here, add some UMAP plots to show the doublets dstribution is desired 
saveRDS(E45.singlet, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E45.singlet.rds")
saveRDS(E55.singlet, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E55.singlet.rds")
saveRDS(E65.singlet, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E65.singlet.rds")
saveRDS(E75.singlet, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E75.singlet.rds")

########################################### end doublefinder  ################################################
################################ https://www.jianshu.com/p/b1947c4156ad   ####################################
#########################   https://github.com/chris-mcginnis-ucsf/DoubletFinder #############################

###########################################  harmony ##########################################
E45.singlet<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E45.singlet.rds")
E55.singlet<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E55.singlet.rds")
E65.singlet<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E65.singlet.rds")
E75.singlet<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\E75.singlet.rds")

library(harmony)

seurat_list <- list(E45.singlet,E55.singlet,E65.singlet,E75.singlet)

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
head(seurat_obj@meta.data)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
######  firstly, evaluating batch effect 
options(repr.plot.height = 5, repr.plot.width = 18)
p1 <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = "time")
p2 <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = "time", pt.size = .1)
plot_grid(p1,p2)

seurat_obj=RunHarmony(seurat_obj,"time", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(seurat_obj, 'harmony')
harmony_embeddings[1:5, 1:5]


options(repr.plot.height = 5, repr.plot.width = 18)
p1 <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = "time")
p2 <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = "time", pt.size = .1)
plot_grid(p1,p2)

# 16. 18. 15.  
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:23)    ## 18
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20,reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

options(repr.plot.width=24,repr.plot.height=4)

DimPlot(seurat_obj,reduction = "umap",label = T)

FeaturePlot(seurat_obj, features = c("PCNA"))
head(seurat_obj@meta.data)


my_color2=c("#C987BB", "#5BBC5E", "#9851A1", "#83BA40", 
            "#8C6DAF", "#ADFF2F", "#5767AE", "#E68824",
            "#6495ED", "#D5AC3C", "#7859A3", "#228B22",
            "#C84F9E", "#4DC1B4", "#EE82EE", "#40E0D0","#2CA029","#FD7F0D")
my_color3=c("#1F8A42", "#832688", "#F47D2C", "#FEE800", 
            "#6E4C9F", "#0C737B", "#282E6C", "#91D5E5",
            "#89C75F", "#E4989C", "#D44C26", "#D7A764",
            "#DB2228", "#87A0D1", "#C06DAD", "#3DBBA7","#2CA029","#FD7F0D")

head(seurat_obj@meta.data)
DimPlot(seurat_obj,reduction = "umap",label = T,cols = my_color3)
DimPlot(seurat_obj,reduction = "umap",label = F,cols = my_color3,group.by = "seurat_clusters")
DimPlot(seurat_obj,reduction = "umap", group.by = "group")


#saveRDS(seurat_obj, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\Integrated_seuratobj.rds")
#seurat_obj<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\Integrated_seuratobj.rds")

FeaturePlot(seurat_obj, features = "ACTA2")
seurat_obj

FeaturePlot(seurat_obj, features = "WNT6", cols = rev(brewer.pal(n = 11, name = "viridis")))


saveRDS(seurat_obj,file = "E:\\Lab of Germ Cell Bio\\10xSpatial_transcriptome\\9.Results\\integrted.rds")
seurat_obj<-readRDS(file = "E:\\Lab of Germ Cell Bio\\10xSpatial_transcriptome\\9.Results\\integrted.rds")
### cell cluster identification 
## markers used 
## Germ cell: Dazl, Stra8, Sycp3             cluster 12
## BPG : Wnt6                                cluster 2,4,14
## EPG:    Krt19 , lhx9                      cluster  7
## Mesenchymal:    Nr2f2, Col1a1             cluster 0,1,3,5,9,10,11,13
## Perivascular:      MMRN1                     cluster 16 
## endothelial     Kdr  Peam1                cluster  8 
## NK/Tcells       CD52                      cluster  15
## Macrophages/Blood related CD74,TYROBP     cluster  6

seurat_obj@meta.data$celltype <- as.numeric(as.character(seurat_obj@meta.data$seurat_clusters))
## assign cell annotation
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==0|celltype==1|celltype==3|celltype==5|celltype==9|celltype==10|celltype==13] <- c("Interstitial"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==11] <- c("Smooth_muscle"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==12] <- c("Germ"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==2|celltype==4|celltype==14] <- c("BPG"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==7] <- c("EPG"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==16] <- c("Perivascular"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==8] <- c("Endothelial"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==15] <- c("T_cell"))
seurat_obj@meta.data <- within(seurat_obj@meta.data, celltype[celltype==6] <- c("Macrophagy"))

DimPlot(seurat_obj,reduction = "umap", group.by = "celltype")
table(seurat_obj@meta.data$celltype)

GR.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GR.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(GR.markers, file = "E:\\Lab of Germ Cell Bio\\10xGenomics_species\\scRNA_pig\\cluster.markers.csv")


## split UMAP plot
levels(seurat_obj)<- c("Germ","BPG","EPG","Interstitial",
                      "Smooth_muscle","Endothelial","Perivascular","T_cell","Macrophagy")
colors <- c("#DB2228","#C06DAD", "#91D5E5","#832688","#DB2228","#89C75F","#2CA029","#3DBBA7","#282E6C")
DimPlot(seurat_obj,reduction = "umap", split.by = "group",cols = colors)
