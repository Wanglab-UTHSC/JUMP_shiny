# ui-covariance.R

fluidPage(
  includeCSS("www/style/theme.css"),
  column(
  3,
  box(title = tagList(icon("tags"), "Select covariate(s) column"),
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