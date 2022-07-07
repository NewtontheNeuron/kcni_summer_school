# This is a script to practice integration with the kcni datasets

library(Seurat)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
human <- read.csv("AIBS_human_counts_mini.csv", row.names = 1)
human_meta <- read.csv("AIBS_human_meta_mini.csv")
mouse <- read.csv("AIBS_mouse_counts_mini.csv", row.names = 1)
mouse_meta <- read.csv("AIBS_mouse_meta_mini.csv")
# Importing the homologous data
homolo <- read.csv("Homologous_genes_mouhum_filtered.csv", row.names = 1)

# Explore and examine the mouse and human data
# first the meta data

head(human)
head(mouse)
head(human_meta)
head(mouse_meta)
head(homolo)

levels(human_meta)
levels(mouse_meta)

# Transpose the matrix
mouse <- t(mouse)
mouse[1:5, 1:5]

# Now merge and coalese
mouse_h <- merge(mouse, homolo[,1:2], by.x = "row.names", by.y = "MGI.symbol", all.x = T, all.y = F)

mouse_h <- mouse_h %>% mutate(HGNC.symbol = coalesce(HGNC.symbol, Row.names))
row.names(mouse_h) <- mouse_h$HGNC.symbol
mouse_h <- mouse_h[, -1]
mouse_h <- mouse_h[, -2870]
mouse <- mouse_h
rm(mouse_h)

# Create the seurat objects
row.names(human_meta) <- human_meta$sample_name
row.names(mouse_meta) <- mouse_meta$sample_name

shuman <- CreateSeuratObject(t(human), meta.data = human_meta)
smouse <- CreateSeuratObject(mouse, meta.data = mouse_meta)

# Create the list and rm the redundant
data.list <- list(smouse, shuman)
rm(human, mouse, human_meta, mouse_meta, shuman, smouse)
gc()

# Save this point so that you do not have start from scratch
saveRDS(data.list, "Data/listofseurats.rds")
data.list <- readRDS("Data/listofseurats.rds")

# Normalize and find variable features for each of the seurat objects separate
data.list <- lapply(data.list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10^6)
  x <- x[homolo$HGNC.symbol,]
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are variable across both seurat objects
features_b <- SelectIntegrationFeatures(object.list = data.list)

anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features_b)
str(anchors)
rm(features_b)

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

# Then follow the regular pipeline.
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 50)
ElbowPlot(combined, ndims = 50)

# Continuing with finding neighbors and clusters
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:35)
combined <- FindClusters(combined, resolution = 0.5)
table(combined@meta.data$seurat_clusters)

# Do the UMAP
combined <- RunUMAP(combined, reduction = "pca", dims = 1:35)
DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap", group.by = "species")

plot2 <- DimPlot(combined, reduction = "umap", group.by = "subclass_label")
plot1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters")
plot1 + plot2

# Save the integration dataset
rm(anchors, data.list, homolo)
gc()
saveRDS(combined, "Data/mouse_human_practice.rds")
combined <- readRDS("Data/mouse_human_practice.rds")

# Data science on clusters
unique(combined$seurat_clusters)
# Creating a mouse and human subclass and species label
combined@meta.data$species_subclass <- paste0(combined@meta.data$species, "_", combined@meta.data$subclass_label)
unique(combined$species_subclass)

as.data.frame.matrix(table(combined$species_subclass, combined$seurat_clusters))
FeaturePlot(combined, features = c("rna_GRIN1", "rna_GRIN2B", "rna_ALDH1L1", "GRIN3A"))

DimPlot(combined, reduction = "umap",
        group.by = "subclass_label",
        cells.highlight = list(combined@meta.data %>%
          filter(subclass_label %in% c("Astro", "Astrocyte")) %>%
          row.names(.),
        combined@meta.data %>%
          filter(subclass_label %in% c("Sst", "SST")) %>%
          row.names(.)),
        cols.highlight = c("lightblue", "goldenrod1"))

# DE between species
# TODO: Use the parallelization
# For performing differential expression after integration, we switch back to the orignal data
DefaultAssay(combined) <- "RNA"
plan("multiprocess", workers = 4) # changes from sequential to parallel

# DE
Idents(combined) <- "seurat_clusters"
cluster_10_markers <- FindConservedMarkers(combined, ident.1 = 10, grouping.var = "species", verbose = F)
head(cluster_10_markers, n = 20)
saveRDS(cluster_10_markers, "Data/astro_human_mouse_markers.rds")
# Limitation - it becomes unitless with the integration so that we cannot just
# look at the log fold change. A three fold change is more vague and more separate

# Differential Expression
DefaultAssay(combined) <- "integrated"

combined$cluster_species <- paste(Idents(combined), combined$species, sep = "_")
Idents(combined) <- "cluster_species"
cluster_10_species_de <- FindMarkers(combined, ident.1 = "10_human", ident.2 = "10_mouse")
head(cluster_10_species_de, n = 25)
saveRDS(cluster_10_species_de, "Data/astro_human_mouse_de.rds")

# Table visualization
library(gt)

head(cluster_10_markers) %>%
  rownames_to_column("genes") %>%
  gt() %>%
  tab_header(
    title = "Mouse and Human Astrocytes"
  ) %>%
  fmt_number(column = -genes, decimals = 2, suffixing = T)
