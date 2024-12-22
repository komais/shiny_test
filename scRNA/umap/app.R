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

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bslib::bs_theme(preset = "bootstrap"),
  # Application title
  titlePanel("单细胞小提琴图可视化"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      datasetInput("rds", file_path, file_name),
      tags$hr(),
      # actionButton("submit", "提交"),
      selectInput("idents", "column", choices = c()),
      colorInput("color")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotlyOutput("distPlot", width = 800, height = 800),
      download_input("download"),
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  read_rds <- datasetServer("rds")

  observeEvent(read_rds(), {
    pbmc <- read_rds()
    columns <- colnames(pbmc@meta.data)
    # 更新下拉菜单的选项
    updateSelectInput(inputId = "idents", choices = columns)
  })

  my_palette <- colorServer("color", read_rds)

  output$distPlot <- renderPlotly({
    req(input$idents, read_rds(), my_palette())
    pbmc <- read_rds()

    idents_length <- length(unique(unlist(pbmc[[input$idents]])))
    if (idents_length <= length(my_palette())) {
      use_palette <- my_palette()[1:idents_length]
    } else {
      use_palette <- my_palette()[1:idents_length %% length(my_palette()) + 1]
    }
    # draw the histogram with the specified number of bins
    p <- DimPlot(
      pbmc,
      reduction = "umap",
      group.by = input$idents,
      cols = use_palette
    )
    ggplotly(p)
  })
  download_serverr("download")
}

# Run the application
shinyApp(ui = ui, server = server)
