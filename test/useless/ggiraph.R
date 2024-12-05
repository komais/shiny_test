library(shiny)
library(ggplot2)
library(ggiraph)
library(dplyr)

# 示例数据
data <- mtcars
data$carname <- rownames(data)

# 定义 UI
ui <- fluidPage(

  # 应用标题
  titlePanel("ggiraph 动态图表示例"),

  # 侧边栏布局
  sidebarLayout(

    # 侧边栏面板
    sidebarPanel(
      selectInput("variable", "选择 X 轴变量:", choices = names(data)[-ncol(data)]),
      sliderInput("pointSize", "点大小:", min = 1, max = 10, value = 3)
    ),

    # 主面板
    mainPanel(
      girafeOutput("distPlot")
    )
  )
)

# 定义服务器逻辑
server <- function(input, output) {

  output$distPlot <- renderGirafe({

    # 根据用户输入过滤数据
    filtered_data <- data

    # 创建 ggplot 图表
    gg_point <- ggplot(data = filtered_data,
                       aes(x = !!sym(input$variable), y = mpg,
                           tooltip = carname, data_id = carname)) +
      geom_point_interactive(aes(size = !!sym(input$variable)), color = "blue", alpha = 0.6) +
      labs(title = paste("MPG vs.", input$variable),
           x = input$variable, y = "MPG") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) # 标题居中


    # 将 ggplot 图表转换为 ggiraph 图表
    girafe(ggobj = gg_point,
           options = list(opts_hover(css = "fill:red;")))
  })
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)
