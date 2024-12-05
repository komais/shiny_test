library(shiny)
library(plotly)

ui <- fluidPage(
  titlePanel("动态 Plotly 图表"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "数据点数量:", min = 1, max = 100, value = 20)
    ),
    mainPanel(
      plotlyOutput("plot")
    )
  )
)

server <- function(input, output) {

  rv <- reactiveValues(x = NULL, y = NULL)

  observeEvent(input$n, {
    rv$x <- 1:input$n
    rv$y <- rnorm(input$n)
  })

  output$plot <- renderPlotly({
    plot_ly(x = rv$x, y = rv$y, type = 'scatter', mode = 'markers')
  })
}

shinyApp(ui = ui, server = server)
