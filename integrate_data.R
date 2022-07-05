library(Seurat)
library(SeuratData)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../../../Neuroscience - CH BSc/NEUR 4908 - F20 W21/RNA Seq/RNAseq (final)/Scripts/")
source("loadLibraries.R")
source("Pre_analysis_functions.R")
mousedata <- load_data("../../clean_neuron_object.RDS")
humandata <- load_data("../../Data/human_ariel_data/top_level_new_annotation.rda")

# Creating the list of Seurat Objects
data.list <- list(mousedata, humandata)

# Now we will follow the introductory integration https://satijalab.org/seurat/articles/integration_introduction.html

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#shuman@assays$RNA@data <- shuman@assays$RNA@counts
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)

anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)
