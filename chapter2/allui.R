library(shiny)
animals <- c("dog", "cat", "mouse", "bird", "other", "I hate animals")
select <- list("1"=c("dog", "cat") , "2"=c("bird", "other"))

ui <- fluidPage(
  textInput("name", "What's your name?", placeholder = "your name"),
  sliderInput("dateto", "When should we deliver", value = as.Date("2020-9-17"), 
              min =as.Date("2020-09-16"), max = as.Date("2020-09-22"),timeFormat="%F"),
  sliderInput("num21", "Number two", value = 50, min = 0, max = 100  , step = 5, animate = T),
  selectInput("ani", "What's your favourite animal?", select),
  
  
  
  passwordInput("password", "What's your password?"),
  textAreaInput("story", "Tell me about yourself", rows = 3),
  numericInput("num", "Number one", value = 0, min = 0, max = 100),
  sliderInput("num2", "Number two", value = 50, min = 0, max = 100),
  sliderInput("rng", "Range", value = c(10, 20), min = 0, max = 100),
  dateInput("dob", "When were you born?"),
  dateRangeInput("holiday", "When do you want to go on vacation next?"),
  selectInput("state", "What's your favourite state?", state.name),
  radioButtons("animal", "What's your favourite animal?", animals),
  radioButtons("rb", "Choose one:",
               choiceNames = list(
                 icon("angry"),
                 icon("smile"),
                 icon("sad-tear")
               ),
               choiceValues = list("angry", "happy", "sad")
  ),
  selectInput(
    "state", "What's your favourite state?", state.name,
    multiple = TRUE
  ),
  fluidRow(
    actionButton("click", "Click me!", class = "btn-danger"),
    actionButton("drink", "Drink me!", class = "btn-lg btn-success")
  ),
  fluidRow(
    actionButton("eat", "Eat me!", class = "btn-block")
  ),
  textOutput("text"),
  verbatimTextOutput("code")
)

server <- function(input, output) {
  output$text <- renderText({ 
    "Hello friend!" 
  })
  output$code <- renderPrint({ 
    summary(1:10) 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)