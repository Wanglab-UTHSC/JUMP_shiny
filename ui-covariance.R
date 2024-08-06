# ui-covariance.R

fluidPage(
  includeCSS("www/style/theme.css"),
  column(
  3,
  box(
    title = tagList(icon("tags"), "Covariates Information"),
    solidHeader = TRUE,
    status = "primary",
    width = NULL,
    textAreaInput(
      "CovInput",
      "Input your covariate group info",
      rows = 6,
      placeholder = paste(
        "Please input covariate information at here. Here is an example format:",
        "-----",
        "sample,sex,age",
        "sig126,F,3",
        "sig127N,M,9",
        "sig127C,M,12",
        "sig128N,F,15",
        "sig128C,F,18",
        "sig129N,M,21",
        "sig129C,F,24",
        "sig130N,M,27",
        sep = '\n'
      )
    ),
    do.call(actionBttn, c(
      list(
        inputId = "confirmedCovList",
        label = "Assign covariates",
        icon = icon("play")),
      actionBttnParams
    )
    ),
    footer = helpText("Please upload the covariance group as example format.")),
  box(title = tagList(icon("tags"), "Select Covariate"),
      solidHeader = TRUE,
      status = "primary",
      width = NULL,
      uiOutput("covVars"),
      do.call(actionBttn, c(
        list(
          inputId = "runCorrection",
          label = "Run Correction",
          icon = icon("play")),
        actionBttnParams
      )
      )
  )
  
),
column(
  9,
  box(
    title = tagList(icon("table"), "Covariates Correction Result Table"),
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    uiOutput("covResultTable")
  ),  
  tabBox(
    title = "",
    width = NULL,
    tabPanel(uiOutput("cor_densityPlot"),title = "Density Plot")
  )
))