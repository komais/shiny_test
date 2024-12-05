library(shiny)
library(reactable)
library(palmerpenguins)

ui <- fluidPage(
  titlePanel("企鹅数据"),
  reactableOutput("penguin_table")
)

server <- function(input, output) {
  output$penguin_table <- renderReactable({
    reactable(
      penguins,
      filterable = TRUE,  # 可过滤
      searchable = TRUE,  # 可搜索
      bordered = TRUE,    # 边框
      striped = TRUE,     # 交替行颜色
      highlight = TRUE,   # 鼠标悬停高亮
      defaultPageSize = 10, # 默认每页显示10行
      columns = list(
        species = colDef(name = "种类"),
        island = colDef(name = "岛屿"),
        bill_length_mm = colDef(name = "喙长 (mm)"),
        bill_depth_mm = colDef(name = "喙宽 (mm)"),
        flipper_length_mm = colDef(name = "鳍长 (mm)"),
        body_mass_g = colDef(name = "体重 (g)"),
        sex = colDef(name = "性别"),
        year = colDef(name = "年份")
      )
    )
  })
}

shinyApp(ui = ui, server = server)
