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