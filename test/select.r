library(shiny)
library(shinyFiles)

ui <- fluidPage(
  shinyFilesButton("file", "选择文件", "请选择一个文件", multiple = FALSE),
  verbatimTextOutput("filepath")
)
server <- function(input, output, session) {
  volumes <- getVolumes() # 获取可用的文件系统卷

  observe({
    shinyFileChoose(input, "file", roots = volumes, session = session)
    
    if (!is.null(input$file)) {
      file_selected <- parseFilePaths(volumes, input$file)
      output$filepath <- renderPrint(file_selected$datapath)
    }
  })
}
shinyApp(ui = ui, server = server)