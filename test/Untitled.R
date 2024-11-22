library(shiny)
# Define UI
ui <- fluidPage(
  # Sidebar with input controls
  sidebarPanel(
    selectInput("groups2", "Select a group variable", choices = c("Group1", "Group2", "Group3")),
    uiOutput("compSelect"),
    uiOutput("group3")
  ),
  # Main panel with output
  mainPanel(
    textOutput("selectedGroups")
  )
)
# Define server logic
server <- function(input, output) {
  # Function to generate example data
  metaData1 <- reactive({
    data.frame(Group1 = sample(c("A", "B", "C"), 20, replace = TRUE),
               Group2 = sample(c("D", "E", "F"), 20, replace = TRUE),
               Group3 = sample(c("G", "H", "I"), 20, replace = TRUE))
  })
  # Render compSelect
  output$compSelect = renderUI({
    group <- input$groups2
    radioButtons(
      inputId = "compareSelection",
      label = "Select pairwise or multiple group comparison",
      choices = c("pairwise" = "pair", "multiple comparison" = "multicomp"),
      selected = "pair"
    )
  })
  # Render group3
  output$group3 <- renderUI({
    output = tagList()
    group <- input$groups2
    df = metaData1()
    if (is.null(colnames(df)) || is.null(group)) {
      vars = NULL
    } else {
      idx = grep(group, colnames(df))
      vars = unique(df[,idx])
    }
    if(input$compareSelection == "pair"){
      output[[1]] = selectInput(
        inputId = "group1select",
        label = "Select first group for comparison",
        choices = vars,
        selected = vars[1]
      )
      output[[2]] = selectInput(
        inputId = "group2select",
        label = "Select second group for comparison",
        choices = vars,
        selected = vars[2]
      )
    }else{
      output[[1]] = selectInput(
        inputId = "group1select",
        label = "Select first group for comparison",
        choices = vars,
        selected = vars[1]
      )
      output[[2]] = selectInput(
        inputId = "group2select",
        label = "Select second group for comparison",
        choices = vars,
        selected = vars[2]
      )
      output[[3]] = selectInput(
        inputId = "group3select",
        label = "Select third group for comparison",
        choices = vars,
        selected = vars[3]
      )
    }
    output
  })
  # Render selectedGroups
  output$selectedGroups <- renderText({
    paste0("You have selected ",
           input$group1select,
           " and ",
           input$group2select,
           if (input$compareSelection == "multicomp") paste0(" and ", input$group3select),
           " for comparison.")
  })
}
# Run the app
shinyApp(ui = ui, server = server) 