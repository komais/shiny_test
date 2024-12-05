library(shiny)
library(Seurat)
library(ggplot2)
library(glue)
library(ggiraph)


#pbmc <- readRDS("pbmc.rds") 

#columns <- colnames(pbmc@meta.data)
datadir <- "E:\\4.bin\\2.R\\4.shiny\\scRNA"
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
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot" , brush="plot_brush" , width="400px"),
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
        output$distPlot <- renderPlot({
        # draw the histogram with the specified number of bins
            DimPlot(pbmc, reduction = "umap", group.by = input$idents)

        } , res=96)

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
                ggsave(file, plot = last_plot(), device = picture_type())
            }
        )
    })
 
}

# Run the application 
shinyApp(ui = ui, server = server)
