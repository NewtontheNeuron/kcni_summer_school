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
