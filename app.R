library(shiny)
library(caret)  # Load the caret package for model handling
library(shinythemes)

# Load the trained RandomForest model (RF8 - RNA-seq)
model_rf8 <- readRDS("AZA_VEN_RandomForest8.rds")  # This is a caret model, not ranger

# Define the User Interface (UI)
ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Apply a modern theme
  titlePanel("AML CR/CRi Probability Prediction Tool (RF8 - RNA-seq)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("Input RNA-seq Data (TPM)", style = "color: #00334d;"),
      numericInput("pdzk1ip1", "PDZK1IP1 (TPM):", value = 0),
      numericInput("afap1", "AFAP1 (TPM):", value = 0),
      numericInput("maneal", "MANEAL (TPM):", value = 0),
      numericInput("sash1", "SASH1 (TPM):", value = 0),
      numericInput("rogdi", "ROGDI (TPM):", value = 0),
      numericInput("cebpe", "CEBPE (TPM):", value = 0),
      numericInput("slc9a3", "SLC9A3 (TPM):", value = 0),
      numericInput("pink1", "PINK1 (TPM):", value = 0),
      
      actionButton("predict_rf8", "Predict CR/CRi Probability", class = "btn-primary", style = "width: 100%; font-weight: bold;")
    ),
    
    mainPanel(
      tags$h2(style = "color: #004c6d; font-weight: bold; text-align: center;", 
              "AML CR/CRi Probability Prediction (RF8 - RNA-seq)"),
      tags$h4("Prediction Results:", style = "color: #004c6d; font-weight: bold; text-align: center;"),
      tags$div(style = "background-color: #f9f9f9; padding: 15px; border-radius: 10px; text-align: center;",
               htmlOutput("result_rf8")),
      tags$p(
        "Disclaimer: This tool is provided for informational purposes only and should NOT be considered as medical advice.",
        style = "font-size: 16px; color: red; font-weight: bold; text-align: center; margin-top: 20px;"
      )
    )
  )
)

# Define the server logic
server <- function(input, output) {
  # Prediction for RF8 (RNA-seq)
  prediction_rf8 <- eventReactive(input$predict_rf8, {
    new_data_rf8 <- data.frame(
      PDZK1IP1 = log2(input$pdzk1ip1 + 1),
      AFAP1 = log2(input$afap1 + 1),
      MANEAL = log2(input$maneal + 1),
      SASH1 = log2(input$sash1 + 1),
      ROGDI = log2(input$rogdi + 1),
      CEBPE = log2(input$cebpe + 1),
      SLC9A3 = log2(input$slc9a3 + 1),
      PINK1 = log2(input$pink1 + 1)
    )
    
    # Use caret's predict() function to get probabilities
    predictions_rf8 <- predict(model_rf8, newdata = new_data_rf8, type = "prob")
    
    # Extract the probability of achieving CR/CRi (class "CR")
    prob_rf8 <- predictions_rf8[, "CR"]  
    
    timestamp_rf8 <- Sys.time()
    
    return(list(prob = prob_rf8, time = timestamp_rf8, data = new_data_rf8))
  })
  
  output$result_rf8 <- renderUI({
    res_rf8 <- prediction_rf8()
    if (!is.null(res_rf8$prob)) {
      result_text <- paste0(
        "The predicted probability of achieving CR/CRi is: ",
        "<span style='color: red; font-weight: bold;'>", 
        round(res_rf8$prob * 100, 2), "%</span>."
      )
      explanation <- paste0(
        "<br><br><strong>Explanation:</strong><br>",
        "This result predicts how likely the patient is to respond to AZA+VEN therapy based on RNA-seq data. ",
        "The likelihood means that out of 100 patients with similar RNA-seq profiles, approximately ", 
        "<span style='font-weight: bold;'>", round(res_rf8$prob * 100, 2), "</span> may benefit from this therapy.<br>",
        "However, this is a rough estimate. Patient outcomes may vary, and a healthcare provider can provide a more precise assessment."
      )
      
      input_summary <- paste0(
        "<br><br><strong>RNA-seq Gene Expression Data (TPM):</strong><br>",
        "PDZK1IP1: ", input$pdzk1ip1, "<br>",
        "AFAP1: ", input$afap1, "<br>",
        "MANEAL: ", input$maneal, "<br>",
        "SASH1: ", input$sash1, "<br>",
        "ROGDI: ", input$rogdi, "<br>",
        "CEBPE: ", input$cebpe, "<br>",
        "SLC9A3: ", input$slc9a3, "<br>",
        "PINK1: ", input$pink1, "<br>",
        "<br><strong>Prediction Time:</strong> ", format(res_rf8$time, "%Y-%m-%d %H:%M:%S")
      )
      
      HTML(paste(result_text, explanation, input_summary))
    } else {
      HTML("<span style='color: red;'>Please input valid data and press 'Predict'.</span>")
    }
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
