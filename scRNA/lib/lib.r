library(shiny)
library(Seurat)
library(ggplot2)
library(glue)
library(plotly)
library(RColorBrewer) # 或其他你re选择的调色板包
library(viridis)
library(purrr)



datasetInput <- function(id, file_path, file_name) {
  file_name_list <- setNames(file_path, file_name)
  selectInput(NS(id, "rds"), "选择rds", choices = file_name_list)
}

datasetServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactive(readRDS(input$rds))
  })
}

download_input <- function(id) {
  taglist <- list(
    selectInput(NS(id, "pic_type"), "Plot Type", choices = c("png", "pdf")),
    downloadButton(NS(id, "download"), "Download Plot")
  )
}

download_serverr <- function(id) {
  moduleServer(id, function(input, output, session) {
    picture_type <- reactive({
      if (input$pic_type == "png") {
        "png"
      } else {
        "pdf"
      }
    })
    output$download <- downloadHandler(
      filename = function() {
        glue("plot.{picture_type()}")
      },
      content = function(file) {
        ggsave(
          file,
          device = picture_type(),
        )
      }
    )
  })
}


colorInput <- function(id) {
  color_sets <- c("Set1", "Set2", "Set3", "Paired", "viridis", "npg", "自定义")
  selectInput(
    NS(id, "palette"), "选择调色板",
    choices = color_sets
  )
}

colorServer <- function(id, read_rds) {
  moduleServer(id, function(input, output, session) {
    #    eventReactive(input$palette, {
    reactive({
      message(input$palette)
      if (input$palette == "自定义") {
        tmp_color <- c(
          "#E41A1C",
          "#377EB8",
          "#3da33a",
          "#984EA3",
          "#FF7F00",
          "#FFFF33",
          "#A65628",
          "#F781BF",
          "#999999",
          "#4C78B8"
        )
      } else if (input$palette == "viridis") {
        tmp_color <- viridis(n = 10)
      } else if (input$palette == "npg") {
        tmp_color <- ggsci::pal_npg()(10)
      } else {
        tmp_color <- brewer.pal(n = 8, name = input$palette)
      }
      message(tmp_color)
      # browser()
      tmp_color
      # my_palette
    })
  })
}
