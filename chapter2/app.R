library(ggplot2)
library(plotly)
library(shiny)

data <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  group = factor(sample(LETTERS[1:3], 100, replace = TRUE))
)

plot <- ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point() +
  labs(title = "Dim Plot in Shiny", x = "X Variable", y = "Y Variable")

plotly_plot <- ggplotly(plot)

ui <- fluidPage(
  plotlyOutput("dimPlot")
)

server <- function(input, output) {
  output$dimPlot <- renderPlotly({
    plotly_plot
  })
}

shinyApp(ui = ui, server = server)