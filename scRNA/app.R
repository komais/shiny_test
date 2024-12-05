#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(Seurat)
library(ggplot2)
library(glue)

pbmc <- readRDS("pbmc.rds") 
browser()


columns <- colnames(pbmc@meta.data)


# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bslib::bs_theme(preset = "bootstrap"),
  
    # Application title
    titlePanel("单细胞可视化"),

    # Sidebar with a slider input for number of bins 
  
    sidebarLayout(
        sidebarPanel(
            selectInput("idents", "column", choices = columns),
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot" , brush="plot_brush"),
           selectInput("pic_type", "Plot Type", choices = c("png", "pdf")),
           downloadButton("download", "Download Plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  thematic::thematic_rmd()
    
    output$distPlot <- renderPlot({

        # draw the histogram with the specified number of bins
        DimPlot(pbmc, reduction = "umap", group.by = input$idents)

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
            ggsave(file, plot = last_plot(), device = picture_type())
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
