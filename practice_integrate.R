# This is a script to practice integration with the kcni datasets

library(Seurat)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
human <- read.csv("AIBS_human_counts_mini.csv")
human_meta <- read.csv("AIBS_human_meta_mini.csv")
mouse <- read.csv("AIBS_mouse_counts_mini.csv")
mouse_meta <- read.csv("AIBS_mouse_meta_mini.csv")

# Explore and examine the mouse and human data
# first the meta data

head(human)
head(mouse)
head(human_meta)
head(mouse_meta)

levels(human_meta)
levels(mouse_meta)