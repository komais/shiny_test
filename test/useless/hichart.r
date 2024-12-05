library(shiny)
library(highcharter)

# 定义UI布局
ui <- fluidPage(
  titlePanel("使用Highcharter生成交互式图表"),
  
  # 添加一个highchart输出
  highchartOutput("myChart")
)

# 定义服务器逻辑
server <- function(input, output) {
  # 创建一个示例数据框
  sample_data <- data.frame(
    Name = c("Alice", "Bob", "Charlie", "David"),
    Age = c(25, 30, 35, 40),
    City = c("New York", "Los Angeles", "Chicago", "Houston")
  )
  
  # 渲染highchart图表
  output$myChart <- renderHighchart({
    hchart(sample_data, "bar", hcaes(x = Name, y = Age)) %>%
      hc_title(text = "Interactive Bar Chart") %>%
      hc_tooltip(pointFormat = "{series.name}: <b>{point.y}</b><br/>Name: {point.x}")
  })
}

# 运行Shiny应用
shinyApp(ui, server)



