library(shiny)
library(colourpicker)

ui <- fluidPage(
  titlePanel("Shiny 取色器示例"),
  sidebarLayout(
    sidebarPanel(
      colourInput("color", "选择颜色", value = "blue")
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    plot(1:10, col = input$color, pch = 16, cex = 2)
  })
}

shinyApp(ui = ui, server = server)