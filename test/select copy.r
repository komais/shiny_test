library(shiny)

ui <- fluidPage(
  actionButton("chooseFile", "选择文件"),
  verbatimTextOutput("filepath")
)

server <- function(input, output, session) {
  filepath <- reactiveVal()

  observeEvent(input$chooseFile, {
    filepath(file.choose())
  })

  output$filepath <- renderPrint({ filepath() })
}

shinyApp(ui = ui, server = server)