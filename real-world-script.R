#Set up Seurat object
library(dplyr)
library(Seurat)
library(patchwork)

#Load data set
merged_counts.data <- read.table('data/merged_counts.txt',sep='\t',header=T)

#Delete entries with ERCC
merged_counts.data <- merged_counts.data[!grepl("ERCC-", merged_counts.data$gene_name),]


#Delete entries of duplicate gene_name
merged_counts.data <- merged_counts.data[!duplicated(merged_counts.data$gene_name),]

#Let gene_name column = row.names
row.names(merged_counts.data) <- merged_counts.data[,2]

#Delete gene_id and gene_name columns
merged_counts.data = subset(merged_counts.data, select = -c(gene_id, gene_name))

#Create Seurat object with raw data
merged_counts <- CreateSeuratObject(counts = merged_counts.data, min.cells = 3, min.features = 200)

#---PRE-PROCESSING---

#Create new column to store QC stats, using set of genes starting with "mt"
merged_counts[["percent.mt"]] <- PercentageFeatureSet(merged_counts, pattern = "^mt-")

head(merged_counts@meta.data, 10)

#Visualize QC metrics
violin_plot <- VlnPlot(merged_counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
feature_plot_1 <- FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "percent.mt")
feature_plot_2 <- FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

violin_plot + feature_plot_1 + feature_plot_2

#Filter data set
merged_counts <- subset(merged_counts, subset = nFeature_RNA > 200 & percent.mt < 5)

#---DATA NORMALIZATION---

#Normalize feature expression measurements for each cell
merged_counts <- NormalizeData(merged_counts)

#---FEATURE SELECTION---
merged_counts <- FindVariableFeatures(merged_counts, selection.method = "vst", nfeature = 2000)

#Get top 10 most highly variable genes
top10 <- head(VariableFeatures(merged_counts), 10)

#Visualize variable features
variable_features_plot_1 <- VariableFeaturePlot(merged_counts)
variable_features_plot_2 <- LabelPoints(plot = variable_features_plot_1, points = top10, repel = TRUE)
variable_features_plot_1 + variable_features_plot_2

#---DATA SCALING---
all.genes <- rownames(merged_counts)
merged_counts <- ScaleData(merged_counts, features = all.genes)

#---LINEAR DIM REDUCTION---
merged_counts <- RunPCA(merged_counts, features = VariableFeatures(object = merged_counts))
print(merged_counts[["pca"]], dims = 1:5, nfeatures = 5)

#Visualize PCs
DimPlot(merged_counts, reduction = "pca")

#---DIMENSIONALITY---

#Visualize dimensionality of data set
ElbowPlot(merged_counts)

#---CLUSTERING---
merged_counts <- FindNeighbors(merged_counts, dims = 1:10)
merged_counts <- FindClusters(merged_counts, resolution = 0.5)

#See cluster IDs of first 5 cells
head(Idents(merged_counts), 50)

#Visualize PCs after clustering
DimPlot(merged_counts, reduction = "pca")

#---NONLINEAR DIM REDUCTION---
merged_counts <- RunUMAP(merged_counts, dims = 1:10)

DimPlot(merged_counts, reduction = "umap", label = TRUE) + NoLegend()

#---CLUSTER BIOMARKERS--- 

#Get top 50 marker genes in each cluster
cluster0.markers <- FindMarkers(merged_counts, 0)
head(cluster0.markers, 50)

cluster1.markers <- FindMarkers(merged_counts, 1)
head(cluster1.markers, 50)

cluster2.markers <- FindMarkers(merged_counts, 2)
head(cluster2.markers, 50)

cluster3.markers <- FindMarkers(merged_counts, 3)
head(cluster3.markers, 50)

cluster4.markers <- FindMarkers(merged_counts, 4)
head(cluster4.markers, 50)

cluster5.markers <- FindMarkers(merged_counts, 5)
head(cluster5.markers, 50)

cluster6.markers <- FindMarkers(merged_counts, 6)
head(cluster6.markers, 50)

cluster7.markers <- FindMarkers(merged_counts, 7)
head(cluster7.markers, 50)

clusterL.markers <- FindMarkers(merged_counts, ident.1=c(0,1,2,3,4), ident.2=6)
clusterL.markers

cluster1.2.markers <- FindMarkers(merged_counts, ident.1=1, ident.2=2)
cluster1.2.markers

cluster2.1.markers <- FindMarkers(merged_counts, ident.1=2, ident.2=1)
cluster2.1.markers


#---DATA VISUALIZATION---

#Generate feature plots for specific genes
FeaturePlot(merged_counts, "Snap25")
FeaturePlot(merged_counts, "Gad1")
FeaturePlot(merged_counts, "Slc17a7")
FeaturePlot(merged_counts, "Slc17a6")

#Generate feature plots for specific genes found in clusters
#From cluster 1
FeaturePlot(merged_counts, "Car4")
FeaturePlot(merged_counts, "Eps8l2")
FeaturePlot(merged_counts, "Rgs6")
FeaturePlot(merged_counts, "Vamp1")
FeaturePlot(merged_counts, "Gsg1l")

#From cluster 2
FeaturePlot(merged_counts, "Necab1")
FeaturePlot(merged_counts, "Cpne6")
FeaturePlot(merged_counts, "Calb1")
FeaturePlot(merged_counts, "Calb2")
FeaturePlot(merged_counts, "6330403K07Rik")

#From cluster L relative to cluster 6
FeaturePlot(merged_counts, "Ntm")
FeaturePlot(merged_counts, "Pcp4")
FeaturePlot(merged_counts, "Cpne7")
FeaturePlot(merged_counts, "Col25a1")
FeaturePlot(merged_counts, "Cbln1")
FeaturePlot(merged_counts, "Kcnq3")
FeaturePlot(merged_counts, "Stum")

#From cluster 1 relative to cluster 2
FeaturePlot(merged_counts, "Cplx1")
FeaturePlot(merged_counts, "Col27a1")
FeaturePlot(merged_counts, "Dpy19l1")
FeaturePlot(merged_counts, "Thsd7a")
FeaturePlot(merged_counts, "Rab37")
FeaturePlot(merged_counts, "Patj")

#From cluster 2 relative to cluster 1
FeaturePlot(merged_counts, "Pcdh10")
FeaturePlot(merged_counts, "Necab1")
FeaturePlot(merged_counts, "Zcchc12")
FeaturePlot(merged_counts, "Scn3b")
FeaturePlot(merged_counts, "Cpne6")


FeaturePlot(merged_counts, "Sec61a1")
FeaturePlot(merged_counts, "Ramp3")

#Data visualization for cluster 1 and cluster 2 with ggplot2
FeaturePlot(merged_counts, c("Dpy19l1", "Necab1"), cols=c("lightgrey", "darkgreen", "purple"),  blend=TRUE) 

#---GGPLOT2 SUMMARY VISUALIZATIONS---

# Import ggplot library
library("tidyverse")

# Find number of cells in each cluster
cells.cluster <- table(Idents(merged_counts))

# Create subsets of cells by cluster
merged_counts.0 <- subset(merged_counts, idents = 0)
merged_counts.1 <- subset(merged_counts, idents = 1)
merged_counts.2 <- subset(merged_counts, idents = 2)
merged_counts.3 <- subset(merged_counts, idents = 3)
merged_counts.4 <- subset(merged_counts, idents = 4)

# Find number of Dpy19l1 in each cluster
Dpy19l1.cluster.1 <- sum(GetAssayData(object = merged_counts.1, slot="data")["Dpy19l1",]>0)
Dpy19l1.cluster.3 <- sum(GetAssayData(object = merged_counts.3, slot="data")["Dpy19l1",]>0)
Dpy19l1.cluster.0 <- sum(GetAssayData(object = merged_counts.0, slot="data")["Dpy19l1",]>0)
Dpy19l1.cluster.4 <- sum(GetAssayData(object = merged_counts.4, slot="data")["Dpy19l1",]>0)
Dpy19l1.cluster.2 <- sum(GetAssayData(object = merged_counts.2, slot="data")["Dpy19l1",]>0)

# Calculate percentage of cells with Dpy19l1 in each cluster
Dpy19l1.1 <- Dpy19l1.cluster.1/cells.cluster[2] * 100
Dpy19l1.3 <- Dpy19l1.cluster.3/cells.cluster[4] * 100
Dpy19l1.0 <- Dpy19l1.cluster.0/cells.cluster[1] * 100
Dpy19l1.4 <- Dpy19l1.cluster.4/cells.cluster[5] * 100
Dpy19l1.2 <- Dpy19l1.cluster.2/cells.cluster[3] * 100

# Find number of Necab1 in each cluster
Necab1.cluster.2 <- sum(GetAssayData(object = merged_counts.2, slot="data")["Necab1",]>0)
Necab1.cluster.4 <- sum(GetAssayData(object = merged_counts.4, slot="data")["Necab1",]>0)
Necab1.cluster.0 <- sum(GetAssayData(object = merged_counts.0, slot="data")["Necab1",]>0)
Necab1.cluster.3 <- sum(GetAssayData(object = merged_counts.3, slot="data")["Necab1",]>0)
Necab1.cluster.1 <- sum(GetAssayData(object = merged_counts.1, slot="data")["Necab1",]>0)

# Calculate percentage of cells with Necab1 in each cluster
Necab1.2 <- Necab1.cluster.2/cells.cluster[3] * 100
Necab1.4 <- Necab1.cluster.4/cells.cluster[5] * 100
Necab1.0 <- Necab1.cluster.0/cells.cluster[1] * 100
Necab1.3 <- Necab1.cluster.3/cells.cluster[4] * 100
Necab1.1 <- Necab1.cluster.1/cells.cluster[2] * 100

# Create dataframe for ggplot
gene_count.data <- data.frame(
  cluster_id = c(1,3,0,4,2),
  Dpy19l1_count = c(Dpy19l1.1, Dpy19l1.3, Dpy19l1.0, Dpy19l1.4, Dpy19l1.2),
  Necab1_count = c(Necab1.1, Necab1.3, Necab1.0, Necab1.4, Necab1.2)
)

# Lock in dataframe order
gene_count.data$cluster_id <- factor(gene_count.data$cluster_id, levels = gene_count.data$cluster_id)

# Create ggplot
library(forcats)
ggplot(gene_count.data, aes(cluster_id, group = 1)) + 
  geom_point(size = 3, aes(y = Dpy19l1_count, color = "Dpy19l1"))+
  geom_line(size = 1, aes(y = Dpy19l1_count, color = "Dpy19l1")) +
  geom_point(size = 3, aes(y = Necab1_count, color = "Necab1"))+
  geom_line(size = 1, aes(y = Necab1_count, color = "Necab1")) +
  theme_classic() +
  ggtitle("Gene expression across endpoints of cluster \'L\'") +
  xlab("Cluster ID") +
  ylab("Percentage of gene expressed") +
  labs(fill = "sd") +
  scale_color_manual(name = "Gene Markers",
                     labels = c("Dpy19l1", "Necab1"),
                     values = c("purple", "darkgreen"))
  



