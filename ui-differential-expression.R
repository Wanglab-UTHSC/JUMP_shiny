#ui-differential-expression.R

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(column(
  3,
  box(
    title = tagList(icon("cogs"), "Differential Expression Methods"),
    id = "DEselect",
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    selectInput(
      "diffMethod",
      "Differential Expression Method",
      choices = c("LIMMA" = "limma")
    ),
    uiOutput("groups2"),
    uiOutput("compSelect"),
    uiOutput("group3"),
    uiOutput("selectInputContainer"),
    # uiOutput("groups3"),
    # uiOutput("groups4"),
    # selectInput("metric2", label = "Select the measure of significance", choice = list("p-value" = "p-value", "FDR" = "FDR"), selected = 1),
    # numericInput("cutoff2", label = "Significance level", min = 0, max = 1, step = 0.01, value = 0.05),
    # numericInput("logfc2", label = "Log2-fold cutoff", value = 1),
    do.call(actionBttn, c(
      list(
        inputId = "DETestType",
        label = "Run Differental Expression Test",
        icon = icon("play")),
      actionBttnParams
    )
    ),
  )
),
column(
  9,
  box(
    title = tagList(icon("table"), "Differential Expression Result Table"),
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    uiOutput("DEResultTable")
  ), 
  # tabsetPanel(
  #   id = "plots",
  #   tabPanel("Volcano Plot", uiOutput("DE_volcanoPlot")),
  #   tabPanel("Heatmap", uiOutput("DE_dendUI")),
  #   tabPanel("Density Plot", uiOutput("DE_distributionUI"))
  # )
  tabBox(
    title = "",
    width = NULL,
    tabPanel(
      title = tagList(icon("bar-chart"), "Moving SD"),
      uiOutput("DE_distributionUI")
    ),
    tabPanel(
      title = tagList(icon("bar-chart"), "SD within group"),
      uiOutput("DE_OldSDUI")
    ),
    tabPanel(
      title = tagList(icon("bar-chart"), "Volcano Plot"),
      uiOutput("DE_volcanoPlot")
    ),
    tabPanel(title = tagList(icon("sitemap"), "Hierarchical Clustering"),
             uiOutput("DE_dendUI"))


  )
)
))