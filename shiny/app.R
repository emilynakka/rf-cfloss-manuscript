library(shiny)
library(tidyverse)
library(ggplot2)

load("rfcfloss_clean_190116.RData")

cfl <- dplyr::rename(cfl, fit = allo)

ui <- fluidPage(
  
  titlePanel(
    title = h2("Regulatory Focus, Counterfactual Loss, and Risky Decision Making: Model Predictions", align="left"),
    windowTitle = "Regulatory Focus, Counterfactual Loss, and Risky Decision Making: Model Predictions"
  ),
  
  sidebarLayout(position = "right",
                sidebarPanel(
                  selectInput("condition", "Condition:", 
                              choices=c("CF Loss","Control")),
                  sliderInput("prev", "Prevention Pride:",
                              min = 1, max = 5,
                              value = 3, step = 0.1),
                  sliderInput("prom", "Promotion Pride:",
                              min = 1, max = 5,
                              value = 3, step = 0.1)
                ),
                mainPanel(
                  htmlOutput("predstatement"),
                  br(),
                  plotOutput("plot2"),
                  br(),
                  br(),
                  p(span("References", style = "font-weight:bold")),
                  p("Nakkawita, E., Mathmann, F., & Higgins, E. T. (in press). Does your gain define my loss?: Socially-defined counterfactual loss and prevention-focused decision-making. ",
                    span("Personality and Individual Differences.", style = "font-style:italic"))
                )
  )
)

server <- function(input, output){
  
  rfcfl.promctr <- reactive({
    rfcfl.promctr <- lm(fit ~ prev * condition + prom * condition, data = cfl)
  })
  
  pred.data <- reactive({
    pred.data <- data.frame(prev = seq(1, 5, .1), prom = input$prom, condition = input$condition)
  })
  
  rfcfl.promctr.pred <- reactive({
    rfcfl.promctr.pred <- cbind(pred.data(), 
                                predict(rfcfl.promctr(), pred.data(), 
                                        interval = "confidence", 
                                        type = c("response", "terms")))
  })
  
  user.pred.data <- reactive({
    user.pred.data <- data.frame(prev = input$prev, prom = input$prom, condition = input$condition)
  })
  
  user.pred <- reactive({
    user.pred <- cbind(user.pred.data(), 
                       predict(rfcfl.promctr(), user.pred.data(), 
                               interval = "confidence", 
                               type = c("response", "terms")))
  })
  
  output$predstatement <- renderText({
    print(paste0("<div style = \"font-size:18px\">Based on the input values you have selected, the model predicts an asset allocation of <b><font color = \"red\">",round(user.pred()$fit, digits = 2),"%</font color></b> to Bitcoin (versus savings).</div>"))
  })
  
  output$plot2 <- renderPlot({
    ggplot(data = cfl, aes(x = prev, y = fit)) +
      xlim(1, 5) +
      ylim(0, 100) +
      geom_segment(data = user.pred(), aes(x = 1, y = fit, xend = prev, yend = fit), linetype = "dashed") +
      geom_segment(data = user.pred(), aes(x = prev, y = 0, xend = prev, yend = fit), linetype = "dashed") +
      geom_ribbon(data = rfcfl.promctr.pred(), aes(ymin = lwr, ymax = upr), alpha = .3, fill = "gray") +
      geom_line(data = rfcfl.promctr.pred(), color = "black", size = 1) +
      geom_point(data = user.pred(), color = "red", size = 7) +
      labs(title="Prevention Pride and Counterfactual Loss as Predictors of Bitcoin Allocation\n(Controlling for Promotion Pride and the Interaction\nBetween Promotion Pride and Counterfactual Loss)", x="Prevention Pride", y="Bitcoin Allocation (%)", color = "Condition")
  })
  
}

shinyApp(ui, server)