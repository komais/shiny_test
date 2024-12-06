library(shiny)
library(Seurat)
library(ggplot2)
library(glue)
library(plotly)
library(RColorBrewer) # 或其他你选择的调色板包
library(viridis)

args <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=", "", args[grep("--file=", args)]))
print(paste("脚本路径:", script_path))

#script_path="E:\\4.bin\\2.R\\4.shiny\\scRNA"
datadir <- glue("{script_path}/data")
datalist <- c("pbmc.rds", "pbmc3k.rds", "pbmc_singlet.rds")

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bslib::bs_theme(preset = "bootstrap"),
  
    # Application title
    titlePanel("单细胞小提琴图可视化"),

    # Sidebar with a slider input for number of bins 
  
    sidebarLayout(
        sidebarPanel(
            selectInput("file1", "选择rds", choices = datalist),
            tags$hr(),
            #actionButton("submit", "提交"),
            selectInput("idents", "column", choices = c()),
            selectInput("palette", "选择调色板", choices = c("Set1", "Set2", "Set3", "Paired", "viridis", "npg", "自定义")),

        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("distPlot" , width=800, height = 800  ),
           selectInput("pic_type", "Plot Type", choices = c("png", "pdf")),
           downloadButton("download", "Download Plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output , session) {
  thematic::thematic_rmd()

  observe({
    url_params <- parseQueryString(session$clientData$url_search)
    
    # 输出第一个参数的值
      if (is.null(url_params[["param_name"]])) {
        message("No parameter found.")
      } else {
        #browser()
        message(paste("The value of 'param_name' is:", url_params[["param_name"]]))
        updateSelectInput(session, "file1", selected = url_params[["param_name"]])
      }
  })
  
  
  
    read_rds <- reactive({
        req( input$file1)
        # 读取文件
        rds_file <- file.path(datadir, input$file1)
        readRDS(rds_file)
    })



    observeEvent(read_rds(), {
        pbmc <- read_rds()
        columns <- colnames(pbmc@meta.data)
        # 更新下拉菜单的选项
        updateSelectInput(inputId = "idents", choices = columns)
    })

    observeEvent(c(input$idents , input$file1 ), {
        req(input$idents, input$file1)
        pbmc <- read_rds()
        idents_length <- length(unique(unlist(pbmc[[input$idents]])))
        output$distPlot <- renderPlotly({
        if (input$palette == "自定义") {
            my_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
        } else if (input$palette == "viridis") {
            my_palette <- viridis(n = idents_length)
        } else if (input$palette == "npg") {
            if (idents_length <= 10) {
                my_palette <-  ggsci::pal_npg()(idents_length)
            } else {
                my_palette <-  ggsci::pal_npg()(10)[1:idents_length %% 10 + 1]
            }

        } else {
            my_palette <- brewer.pal(n = idents_length, name = input$palette)
        }
        # draw the histogram with the specified number of bins
            p <- DimPlot(pbmc, reduction = "umap", group.by = input$idents , cols = my_palette)

            ggplotly(p)

        })

        picture_type <- reactive({
            if (input$pic_type == "png") {
                "png"
            } else {
                "pdf"
            }
        })
        # Download button
        output$download <- downloadHandler(
            filename = function() {
                glue("plot.{picture_type()}")
            },
            content = function(file) {
                p <- DimPlot(pbmc, reduction = "umap", group.by = input$idents)
                if (picture_type() == "pdf") {
                    ggsave(file, plot = p, device = picture_type() , width = 10, height = 10)
                } else {
                    ggsave(file, plot = p, device = picture_type(), width = 10, height = 10)
                }
            }
        )
    })
 
}

# Run the application 
app <- shinyApp(ui = ui, server = server )


runApp(app , port=4949)
