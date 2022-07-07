library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
human <- read.csv("AIBS_human_counts_mini.csv")
human_meta <- read.csv("AIBS_human_meta_mini.csv")

head(human_meta)
human_meta[1:3, 1:3]

# Looking at only the neuronal cells.
human_meta_filtered <- human_meta %>%
  filter(class_label != "Non-neuronal") %>%
  select(sample_name) %>%
  unlist()

human[human_meta_filtered, ][1:10]

human %>%
  filter(row.names(.) %in% human_meta_filtered)


#### Day 2
library(Seurat)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
human <- read.csv("AIBS_human_counts_mini.csv", row.names = 1)
human_meta <- read.csv("AIBS_human_meta_mini.csv")

#TODO: Seurat matches the column names of the count matrix to the
# row names in the meta data.
row.names(human_meta) <- human_meta$sample_name
row.names(human_meta)[1:5]

human_tpose <- t(human)
human_meta[1:5, 1:5]
human_tpose[1:5, 1:5]
# Make the Seurat object
shuman <- CreateSeuratObject(counts = human_tpose,
                             meta.data = human_meta)
# TODO: If the counts has decimal points it is likely a normalized
# S3 object is a list of lists. You put together lists
# S4 object is more formal.
# Assay is equivalent to summarized experiments.

# Normalize
shuman <- NormalizeData(shuman, normalization.method = "LogNormalize", scale.factor = 10^6)
shuman@assays$RNA@counts[1:10, 1:4]
shuman@assays$RNA@data[1:10, 1:4]

# Find variable features
# TODO: these genes are the highest highs and lowest lows. They are not for very specific cells having 100 and everyone zero.

shuman <- FindVariableFeatures(shuman, selection.method = "vst", nfeatures = 2000)
shuman@assays$RNA@var.features[1:5]
# TODO: "We have gone higher but never lower 8000 4000 maybe." Sonny. Higher
# to get more consistent clustering
# Do do not want to go too low

# Scaling is so that the variance of highly expressed genes do not overpower
# the lower expressed genes.
shuman <- ScaleData(shuman)

# TODO: PCAs should be used with highly correlated variables.
shuman <- RunPCA(shuman) #, features = shuman@assays$RNA@var.features)
ElbowPlot(shuman, ndims = 50)
# You can in another analysis get a flat line.
# not sure what happens if you force it after 
# This part is more for the picking the dimensional reductions we use in downstream
# analyses

# What do you consistently get versus whether any of these numbers are good or bad

# Clusters
shuman <- FindNeighbors(shuman, dims = 1:20) # based on what we found earlier
shuman <- FindClusters(shuman, resolution = 0.5)
# changing the parameters to the Find* functions will drastically change your results.

shumna$seurat_clusters[1:5]

table(shumna$seurat_clusters, shuman$subclass_label)

# Visulaization
shuman <- RunUMAP(shumna, dims = 1:20)
DimPlot(shuman, reduction = "umap")
DimPlot(shuman, reduction = "pca")
# Everyone gets the same plot because it is reducing to 2 dimensions
Idents(shuman) <- "subclass_label"
DimPlot(shuman, reduction = "umap")
DimPlot(shuman, reduction = "umap", group.by = "seurat_clusters")

saveRDS(shuman, "Data/human_practice.rds")


## Day 3

library(Seurat)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
mouse <- read.csv("AIBS_mouse_counts_mini.csv", row.names = 1)
mouse_meta <- read.csv("AIBS_mouse_meta_mini.csv")

row.names(mouse_meta) <- mouse_meta$sample_name

smouse <- CreateSeuratObject(counts = t(mouse), meta.data = mouse_meta)

smouse <- NormalizeData(smouse, normalization.method = "LogNormalize", scale.factor = 10^6)

smouse <- FindVariableFeatures(smouse, selection.method = "vst", nfeatures = 2000)

smouse <- ScaleData(smouse)

smouse <- RunPCA(smouse, npcs = 45)
ElbowPlot(smouse, ndims = 45)

# Do not do the bellow ##########################

smouse <- FindNeighbors(smouse, dims = 13)
smouse <- FindClusters(smouse, resolution = 0.5)
smouse <- FindNeighbors(smouse, dims = 20)
smouse <- FindClusters(smouse, resolution = 0.5)

# Continue with this #############################

smouse <- FindNeighbors(smouse, dims = 30)
smouse <- FindClusters(smouse, resolution = 0.5)

# Do not do the bellow ###########################

resolutions_to_try <- list(0.1, 0.25, 0.5, 1.2, 2, 5, 8)
smouse <- FindNeighbors(smouse, dims = 30)

smouse_list <- lapply(resolutions_to_try, function(x){
  smouse <- FindClusters(smouse, resolution = x)
  return(smouse)
})

# skip this normally
dimplot <- list()
lapply(seq_along(smouse_list), function(x){
  print(resolutions_to_try[[x]])
  #table(smouse_list[[x]]$seurat_clusters, smouse_list[[x]]$subclass_label)
  smouse_list[[x]] <- RunUMAP(smouse_list[[x]], dims = 1:20)
  Idents(smouse_list[[x]]) <- "subclass_label"
  dimplot[[x]] <<- DimPlot(smouse_list[[x]], reduction = "umap")
})

# Continue down here #############################

table(smouse$seurat_clusters, smouse$subclass_label)

smouse <- RunUMAP(smouse, dims = 1:20)
Idents(smouse) <- "subclass_label"
DimPlot(smouse, reduction = "umap")

saveRDS(smouse, "Data/mouse_practice.rds")


# Day 3 part 2
shuman <- readRDS("Data/human_practice.rds")
Idents(shuman) <- "seurat_clusters"
Idents(shuman) <- "donor_sex_label"
# We want to be able to distinguish between a VIP cell vs LAMP5 cell
DimPlot(shuman, reduction = "umap", group.by = "subclass_label", label = T, repel = T)

# We will look between vip and Pax6 neurons.
vip_pax6_markers <- FindMarkers(shuman, ident.1 = 3, ident.2 = 8, logfc.threshold = 0.25, min.pct = 0.125)

vip_pax6_markers %>%
  arrange(desc(avg_log2FC))

cluster_3_vall <- FindMarkers(shuman, ident.1 = 3, logfc.threshold = 0.25,
                              min.pct = 0.50)
cluster_3_vall %>%
  arrange(desc(avg_log2FC))

cluster_4_vall <- FindMarkers(shuman, ident.1 = 4, logfc.threshold = 0.25,
                              min.pct = 0.50) # look into test.use.
cluster_4_vall %>%
  arrange(desc(avg_log2FC))

cluster_4_vall[which(row.names(cluster_4_vall) == "GRIN1"), ]

# 

features <- cluster_4_vall %>%
  arrange(p_val_adj) %>%
  head(n=6) %>%
  row.names()

features

# Violin Plot
VlnPlot(shuman, features = features)

# Feature Plot
FeaturePlot(shuman, features = features)

# Dot Plot
DotPlot(shuman, features = features) + RotatedAxis()

# Heat Map
DoHeatmap(subset(shuman, downsample = 100), features = features, size = 3, slot = "data") + scale_fill_viridis_c(option = "magma")


# Day 4

