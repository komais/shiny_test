library(shiny)
library(rhandsontable)

ui <- fluidPage(
  rHandsontableOutput("hot")
)

server <- function(input, output) {
  output$hot <- renderRHandsontable({
    df <- data.frame(
      col1 = c(1, 2, 3),
      col2 = c("A", "B", "C")
    )
    rhandsontable(df)
  })
}

shinyApp(ui = ui, server = server)