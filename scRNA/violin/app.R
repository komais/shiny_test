library(shiny)
library(Seurat)
library(ggplot2)
library(glue)
library(plotly)
library(RColorBrewer) # 或其他你re选择的调色板包
library(viridis)
library(purrr)
source("../lib/lib.r")
script_path <- "E:\\4.bin\\2.R\\4.shiny\\scRNA\\data"
file_name <- c("pbmc.rds", "pbmc3k.rds", "pbmc_singlet.rds")
file_path <- map(file_name, function(x) {
  glue("{script_path}\\{x}")
})


