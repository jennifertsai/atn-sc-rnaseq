##############################################################################
# Jennifer Tsai
# sc-RNA-seq Analysis of anterior thalamic nuclei
# January-April 2022
##############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse) #ggplot library

##############################################################################

#Load data set
merged_counts.data <-
  read.table('data/merged_counts.txt', sep = '\t', header = T)

#DATA-PREPROCESSING###########################################################

#Delete entries with ERCC
merged_counts.data <-
  merged_counts.data[!grepl("ERCC-", merged_counts.data$gene_name),]


#Delete entries of duplicate gene_name
merged_counts.data <-
  merged_counts.data[!duplicated(merged_counts.data$gene_name),]

#Let gene_name column = row.names
row.names(merged_counts.data) <- merged_counts.data[, 2]

#Delete gene_id and gene_name columns
merged_counts.data = subset(merged_counts.data, select = -c(gene_id, gene_name))

#Create Seurat object with raw data
merged_counts <-
  CreateSeuratObject(counts = merged_counts.data,
                     min.cells = 3,
                     min.features = 200)

#Create new column to store QC stats, using set of genes starting with "mt"
merged_counts[["percent.mt"]] <-
  PercentageFeatureSet(merged_counts, pattern = "^mt-")

head(merged_counts@meta.data, 10)

#Visualize QC metrics
violin_plot <-
  VlnPlot(
    merged_counts,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )
feature_plot_1 <-
  FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "percent.mt")
feature_plot_2 <-
  FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

violin_plot + feature_plot_1 + feature_plot_2

#Filter data set
merged_counts <-
  subset(merged_counts, subset = nFeature_RNA > 1500 &
           percent.mt < 5)

#DATA NORMALIZATION###########################################################

#Normalize feature expression measurements for each cell
merged_counts <- SCTransform(merged_counts)
#normalize with SCT

#FEATURE SELECTION############################################################

merged_counts <-
  FindVariableFeatures(merged_counts,
                       selection.method = "vst",
                       nfeature = 2000)

#Get top 10 most highly variable genes
top10 <- head(VariableFeatures(merged_counts), 10)

#Visualize variable features
variable_features_plot_1 <- VariableFeaturePlot(merged_counts)
variable_features_plot_2 <-
  LabelPoints(plot = variable_features_plot_1,
              points = top10,
              repel = TRUE)
variable_features_plot_1 + variable_features_plot_2

#DATA SCALING#################################################################

all.genes <- rownames(merged_counts)
merged_counts <- ScaleData(merged_counts, features = all.genes)

#LINEAR DIM REDUCTION#########################################################

merged_counts <-
  RunPCA(merged_counts, features = VariableFeatures(object = merged_counts))
print(merged_counts[["pca"]], dims = 1:5, nfeatures = 5)

#Visualize PCs
DimPlot(merged_counts, reduction = "pca")

#DIMENSIONALITY###############################################################

#Visualize dimensionality of data set
ElbowPlot(merged_counts)

#CLUSTERING###################################################################

merged_counts <- FindNeighbors(merged_counts, dims = 1:10)
merged_counts <- FindClusters(merged_counts, resolution = 0.10)

#See cluster IDs of first 5 cells
head(Idents(merged_counts), 50)

#Visualize PCs after clustering
DimPlot(merged_counts, reduction = "pca")

#Visualize QC metrics after clustering
violin_plot <-
  VlnPlot(
    merged_counts,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )
feature_plot_1 <-
  FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "percent.mt")
feature_plot_2 <-
  FeatureScatter(merged_counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

violin_plot + feature_plot_1 + feature_plot_2

#NONLINEAR DIM REDUCTION######################################################

merged_counts <- RunUMAP(merged_counts, dims = 1:10)

umap.plot <-
  DimPlot(merged_counts, reduction = "umap", label = TRUE) + NoLegend()

umap.plot

#CLUSTER BIOMARKERS###########################################################

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

#Gene markers for L shaped clump in UMAP
clusterL.markers <-
  FindMarkers(merged_counts,
              ident.1 = c(0, 1, 2, 3),
              ident.2 = 4)
clusterL.markers

#Gene markers of clusters relative to each another cluster
cluster3.0.markers <-
  FindMarkers(merged_counts, ident.1 = 3, ident.2 = 0)
head(cluster3.0.markers, 50)

cluster0.3.markers <-
  FindMarkers(merged_counts, ident.1 = 0, ident.2 = 3)
head(cluster0.3.markers, 50)

cluster2.3.markers <-
  FindMarkers(merged_counts, ident.1 = 2, ident.2 = 3)
head(cluster2.3.markers, 50)

#DATA VIZ#####################################################################

#Generate feature plots for specific genes
FeaturePlot(merged_counts, "Snap25")
FeaturePlot(merged_counts, "Gad1")
FeaturePlot(merged_counts, "Slc17a7")
FeaturePlot(merged_counts, "Slc17a6")

#Generate feature plots for specific genes found in clusters
#From cluster 4 (new)
FeaturePlot(merged_counts, "C1ql2")
FeaturePlot(merged_counts, "Eepd1")
FeaturePlot(merged_counts, "Gng13")
FeaturePlot(merged_counts, "Slc17a7")
FeaturePlot(merged_counts, "Rspo3")


#From cluster 0 gene markers (new)
FeaturePlot(merged_counts, "Car4")
FeaturePlot(merged_counts, "Rab37")
FeaturePlot(merged_counts, "Dpy19l1")
FeaturePlot(merged_counts, "Kcnab3")
FeaturePlot(merged_counts, "Eps8l2")


FeaturePlot(merged_counts,
            c("Dpy19l1",
              "Col27a1",
              "Rab37",
              "Car4",
              "Eps8l2",
              "Rgs6"))
FeaturePlot(merged_counts,
            c("Necab1",
              "Pcdh10",
              "Zcchc12",
              "Cpne6",
              "Scn3b",
              "Calb1"))


#From cluster 3 gene markers (new)
FeaturePlot(merged_counts, "Necab1")
FeaturePlot(merged_counts, "Cpne6")
FeaturePlot(merged_counts, "Calb1")
FeaturePlot(merged_counts, "Calb2")
FeaturePlot(merged_counts, "Pcdh10")

#From cluster L relative to cluster 4
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


FeaturePlot(merged_counts, "Reep5")
FeaturePlot(merged_counts, "Ramp3")

#Data visualization for cluster 1 and cluster 2 with ggplot2
FeaturePlot(
  merged_counts,
  c("Dpy19l1", "Necab1"),
  blend = TRUE,
  cols= c("red", "blue")
) 

#GGPLOT2 SUMMARY VIZ##########################################################

#Find number of cells in each cluster
cells.cluster <- table(Idents(merged_counts))

#Get percentage of gene in each cluster
gene_percent_in_clusters <- function(gene_name) {
  all_genes.percentages <- c()
  for (i in 0:4) {
    merged_counts.cluster <- subset(merged_counts, idents = i)
    gene.cluster <-
      sum(GetAssayData(object = merged_counts.cluster, slot = "data")[gene_name,] >
            0)
    gene.percent <- gene.cluster / cells.cluster[i + 1] * 100
    all_genes.percentages <- c(all_genes.percentages, gene.percent)
  }
  all_genes.percentages
}

#Percent of cluster 1 gene markers
Dpy19l1_percents <- gene_percent_in_clusters("Dpy19l1")
Col27a1_percents <- gene_percent_in_clusters("Col27a1")
Rab37_percents <- gene_percent_in_clusters("Rab37")
Car4_percents <- gene_percent_in_clusters("Car4")
Eps8l2_percents <- gene_percent_in_clusters("Eps8l2")
Rgs6_percents <- gene_percent_in_clusters("Rgs6")

#Percent of cluster 2 gene markers
Necab1_percents <- gene_percent_in_clusters("Necab1")
Pcdh10_percents <- gene_percent_in_clusters("Pcdh10")
Zcchc12_percents <- gene_percent_in_clusters("Zcchc12")
Cpne6_percents <- gene_percent_in_clusters("Cpne6")
Scn3b_percents <- gene_percent_in_clusters("Scn3b")
Calb1_percents <- gene_percent_in_clusters("Calb1")

#Create matrix of percents of all gene markers
percents_list <- t(matrix(
  c(
    Dpy19l1_percents,
    Col27a1_percents,
    Rab37_percents,
    Car4_percents,
    Eps8l2_percents,
    Rgs6_percents,
    Necab1_percents,
    Pcdh10_percents,
    Zcchc12_percents,
    Cpne6_percents,
    Scn3b_percents,
    Calb1_percents
  ),
  nrow = 5,
  ncol = 12
))

#Create dataframe for 'specific-gene' ggplot
gene_count.data <- data.frame(
  cluster_id = c(1, 3, 0, 4, 2),
  Dpy19l1_count = Dpy19l1_percents[c(2, 4, 1, 5, 3)],
  Necab1_count = Necab1_percents[c(2, 4, 1, 5, 3)]
)

#Lock in dataframe order
gene_count.data$cluster_id <-
  factor(gene_count.data$cluster_id, levels = gene_count.data$cluster_id)

#Create ggplot

ggplot(gene_count.data, aes(cluster_id, group = 1)) +
  geom_point(size = 3, aes(y = Dpy19l1_count, color = "Dpy19l1")) +
  geom_line(size = 1, aes(y = Dpy19l1_count, color = "Dpy19l1")) +
  geom_point(size = 3, aes(y = Necab1_count, color = "Necab1")) +
  geom_line(size = 1, aes(y = Necab1_count, color = "Necab1")) +
  theme_classic() +
  ggtitle("Gene expression across endpoints of cluster \'L\'") +
  xlab("Cluster ID") +
  ylab("Percentage of gene expressed") +
  labs(fill = "sd") +
  scale_color_manual(
    name = "Gene Markers",
    labels = c("Dpy19l1", "Necab1"),
    values = c("purple", "darkgreen")
  )

#Create dataframe for 'all-gene' ggplot
gene_names <-
  c(
    "Dpy19l1",
    "Col27a1",
    "Rab37",
    "Car4",
    "Eps8l2",
    "Rgs6",
    "Necab1",
    "Pcdh10",
    "Zcchc12",
    "Cpne6",
    "Scn3b",
    "Calb1"
  )
cluster_order <- c(1, 3, 0, 4, 2)

#Create dataframe for all-gene box plot
genes.percentages <- data.frame()
for (i in cluster_order) {
  for (j in gene_names) {
    if (which(gene_names == j) >= 1 && which(gene_names == j) <= 6) {
      cluster_id <- 1
    }
    else {
      cluster_id <- 2
    }
    genes.percentages <- rbind(
      genes.percentages,
      data.frame(
        gene = j,
        cluster = i,
        id = cluster_id,
        percent_gene = percents_list[which(gene_names == j), i + 1]
      )
    )
  }
}

#Separate cluster 1 and cluster 2 data
cluster1.data <- subset(genes.percentages, id == 1)
cluster2.data <- subset(genes.percentages, id == 2)

cluster1.data$cluster <- as.factor(cluster1.data$cluster)
cluster2.data$cluster <- as.factor(cluster2.data$cluster)

#Create all-gene box plot for cluster 1
cluster1.plot <-
  (
    ggplot(data = cluster1.data, aes(x = fct_inorder(cluster),
                                     y = percent_gene,))
    + geom_boxplot()
    + geom_dotplot(
      binaxis = 'y',
      dotsize = 0.5,
      stackdir = 'center'
    )
    + geom_line(aes(group = gene, color = gene), size = 1)
    + ggtitle("Top 6 cluster 1 gene markers across endpoints of cluster \'L\'")
    + xlab("Cluster ID")
    + ylab("% of gene expressed")
    + theme_classic()
  )

#Create all-gene box plot for cluster 2
cluster2.plot <-
  (
    ggplot(data = cluster2.data, aes(x = fct_inorder(cluster),
                                     y = percent_gene))
    + geom_boxplot()
    + geom_dotplot(
      binaxis = 'y',
      dotsize = 0.5,
      stackdir = 'center'
    )
    + geom_line(aes(group = gene, color = gene), size = 1)
    + ggtitle("Top 6 cluster 2 gene markers across endpoints of cluster \'L\'")
    + xlab("Cluster ID")
    + ylab("% of gene expressed")
    + theme_classic()
  )

cluster1.plot + cluster2.plot

#CONTINUUM PLOT###############################################################

#Generate interactive UMAP to find coordinates of endpoints
HoverLocator(plot = umap.plot,
             information = FetchData(
               merged_counts,
               vars = c("ident", "UMAP_1", "UMAP_2", "PC_1", "PC_2", "nFeature_RNA")
             ))

#Pick Q8_C08_Plate_502_Q8 and Q7_A09_Plate_503_Q7 coordinates as arbitrary seed values
coord_start <- c(-7.756, -0.694)
coord_end <- c(1.583, 7.709)
coord_table <- rbind(coord_start, coord_end)
continuum <- dist(coord_table, "manhattan")

#Create subset with only cells from cluster L
clusterL.data <-
  subset(
    merged_counts,
    seurat_clusters == 0 |
      seurat_clusters == 1 |
      seurat_clusters == 2 |
      seurat_clusters == 3 |
      seurat_clusters == 4
  )
clusterL_coords <- clusterL.data[["umap"]]@cell.embeddings

#Calculate gene expression percentages of each cell
#Cluster 1
percent.avg1 <-
  (
    PercentageFeatureSet(clusterL.data, pattern = "Car4") +
      PercentageFeatureSet(clusterL.data, pattern = "Col27a1") +
      PercentageFeatureSet(clusterL.data, pattern = "Dpy19l1") +
      PercentageFeatureSet(clusterL.data, pattern = "Eps8l2") +
      PercentageFeatureSet(clusterL.data, pattern = "Rab37") +
      PercentageFeatureSet(clusterL.data, pattern = "Rgs6")
  ) / 6

#Cluster 2
percent.avg2 <-
  (
    PercentageFeatureSet(clusterL.data, pattern = "Calb1") +
      PercentageFeatureSet(clusterL.data, pattern = "Cpne6") +
      PercentageFeatureSet(clusterL.data, pattern = "Necab1") +
      PercentageFeatureSet(clusterL.data, pattern = "Pcdh10") +
      PercentageFeatureSet(clusterL.data, pattern = "Scn3b") +
      PercentageFeatureSet(clusterL.data, pattern = "Zcchc12")
  ) / 6

#Intermediate Clusters
percent.avg0 <-
  (
    PercentageFeatureSet(clusterL.data, pattern = "Reep5") +
      PercentageFeatureSet(clusterL.data, pattern = "Slc25a4") +
      PercentageFeatureSet(clusterL.data, pattern = "Rps3a1") +
      PercentageFeatureSet(clusterL.data, pattern = "Stmn1") +
      PercentageFeatureSet(clusterL.data, pattern = "Tuba1a") +
      PercentageFeatureSet(clusterL.data, pattern = "Reep5")
  ) / 6

FeaturePlot(clusterL.data, "Ldhb")
#Get Manhattan distance from starting coordinate to each cell
cell_distances <- c()
normalized_distances <- c()
for (row in 1:nrow(clusterL_coords)) {
  coord_cell <- c(clusterL_coords[row, 1], clusterL_coords[row, 2])
  coords <- rbind(coord_start, coord_cell)
  start_to_cell <- dist(coords, "manhattan")
  cell_distances <- c(cell_distances, start_to_cell)
  
  #Normalized distances with continuum scale
  if (start_to_cell > continuum) {
    normalized_distance <- 1
  }
  else {
    normalized_distance <- start_to_cell / continuum
  }
  normalized_distances <-
    c(normalized_distances, normalized_distance)
}

#Create ggplot dataframe
continuum.data1 <- data.frame(percent.avg1, normalized_distances)
continuum.data2 <- data.frame(percent.avg2, normalized_distances)
continuum.data0 <- data.frame(percent.avg0, normalized_distances)


continuum.data1$nCount_RNA <- unlist(continuum.data1$nCount_RNA)
continuum.data1$normalized_distances <-
  unlist(continuum.data1$normalized_distances)

continuum.data2$nCount_RNA <- unlist(continuum.data2$nCount_RNA)
continuum.data2$normalized_distances <-
  unlist(continuum.data2$normalized_distances)

continuum.data0$nCount_RNA <- unlist(continuum.data0$nCount_RNA)
continuum.data0$normalized_distances <-
  unlist(continuum.data0$normalized_distances)

#Create continuum plots
continuum.plot1 <-
  (
    ggplot(data = continuum.data1, aes(x = normalized_distances, y = nCount_RNA))
    + geom_point(alpha=0.25)
    + theme_classic()
    + ggtitle("Average cluster 1 gene expressions across cluster \'L\' cells")
    + xlab("Distance along continuum")
    + ylab("% of top 6 gene markers expressed")
    
  )

continuum.plot2 <-
  (
    ggplot(data = continuum.data2, aes(x = normalized_distances, y = nCount_RNA))
    + geom_point(alpha=0.25)
    + theme_classic()
    + ggtitle("Average cluster 2 gene expressions across cluster \'L\' cells")
    + xlab("Distance along continuum")
    + ylab("% of top 6 gene markers expressed")
  )

continuum.plot0 <-
  (
    ggplot(data = continuum.data0, aes(x = normalized_distances, y = nCount_RNA))
    + geom_point(alpha=0.25)
    + theme_classic()
    + ggtitle("Average inter. gene expressions across cluster \'L\' cells")
    + xlab("Distance along continuum")
    + ylab("% of top 6 gene markers expressed")
  )

continuum.plot1 + continuum.plot2 +continuum.plot0
