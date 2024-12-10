library(shiny)
library(glue)


problem <- function(seed_offset = 0) {
  set.seed(as.integer(Sys.time()) + seed_offset) # 使用不同的偏移量来确保不同的种子
  n <- 1 # 每次只生成一个问题
  num1 <- sample(0:100, n)
  op <- sample(c("+", "-"), n)
  num2 <- ifelse(op == "+", sample(0:(100 - num1), 1), sample(0:num1, 1))

  correct_answer <- if (op == "+") {
    num1 + num2
  } else {
    num1 - num2
  }

  return(data.frame(
    num1 = num1, op = op, num2 = num2,
    correct_answer = correct_answer
  ))
}

questions <- lapply(0:9, function(i) problem(i)) # 使用0到9作为seed_offset
questions_df <- do.call(rbind, questions) # 将列表合并成一个数据框


# UI 部分
ui <- fluidPage(
  titlePanel("加减法小测验"),
  lapply(1:10, function(i) {
    numericInput(glue("problem{i}"),
      glue("第{i}题:  {questions_df$num1[i]} {questions_df$op[i]} {questions_df$num2[i]} 等于？"),
      value = ""
    )
  }),
  actionButton("submit", "提交"),
  verbatimTextOutput("result")
)

score <- function(list1, list2) {
  total_num <- length(list1)
  correct_num <- 0
  for (i in 1:total_num) {
    if (list1[[i]] == list2[[i]]) {
      correct_num <- correct_num + 1
    }
  }
  return(glue("正确率：{round(correct_num/total_num*100,2)}分"))
}

# Server 部分
server <- function(input, output, session) {
  # 创建一个reactiveValues对象来存储数值
  values <- reactiveValues(list = list())

  # 观察"collect"按钮的点击事件
  observeEvent(input$submit, {
    # 收集所有numericInput组件的值
    values$list <- lapply(1:10, function(i) {
      input[[glue("problem{i}")]]
    })

    # 打印收集到的值到控制台（可选）
    # print(values$list)
  })

  # 输出收集到的值
  output$result <- renderPrint({
    if (length(values$list) == 0) {
      return("请先点击提交按钮")
    } else {
      return(score(values$list, questions_df$correct_answer))
    }
  })
}



# 运行应用
shinyApp(ui = ui, server = server)
