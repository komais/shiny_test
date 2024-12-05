library(shiny)

# 定义UI
ui <- fluidPage(
  titlePanel("100以内加减法练习"),
  sidebarLayout(
    sidebarPanel(
      actionButton("generate", "生成新题目")
    ),
    mainPanel(
      uiOutput("questions"),
      textInput("answer1", label = "第1题答案:", value = ""),
      textInput("answer2", label = "第2题答案:", value = ""),
      textInput("answer3", label = "第3题答案:", value = ""),
      textInput("answer4", label = "第4题答案:", value = ""),
      textInput("answer5", label = "第5题答案:", value = ""),
      textInput("answer6", label = "第6题答案:", value = ""),
      textInput("answer7", label = "第7题答案:", value = ""),
      textInput("answer8", label = "第8题答案:", value = ""),
      textInput("answer9", label = "第9题答案:", value = ""),
      textInput("answer10", label = "第10题答案:", value = ""),
      actionButton("submit", "提交答案"),
      verbatimTextOutput("score")
    )
  )
)

# 定义服务器逻辑
server <- function(input, output, session) {
  # 创建反应值存储问题和答案
  questions <- reactiveValues(problems = NULL, answers = rep(NA, 10))
  
  observeEvent(input$generate, {
    # 生成10个随机的加减法问题
    problems <- lapply(1:10, function(i) {
      a <- sample(1:100, 1)
      b <- sample(1:100, 1)
      op <- sample(c("+", "-"), 1)
      list(a = a, b = b, op = op, answer = ifelse(op == "+", a + b, a - b))
    })
    questions$problems <- problems
  })
  
  output$questions <- renderUI({
    if (is.null(questions$problems)) return()
    lapply(1:10, function(i) {
      div(paste0(questions$problems[[i]]$a, " ", 
                 questions$problems[[i]]$op, " ", 
                 questions$problems[[i]]$b, " = "),
          br())
    })
  })
  
  observeEvent(input$submit, {
    user_answers <- sapply(1:10, function(i) as.numeric(input[[paste0("answer", i)]]))
    correct_answers <- sapply(questions$problems, `[[`, "answer")
    score <- sum(user_answers == correct_answers, na.rm = TRUE)
    output$score <- renderPrint({
      paste("你的得分是：", score, "/10")
    })
  })
}

# 运行应用
shinyApp(ui = ui, server = server)