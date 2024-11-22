# ui-homepage.R

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(
  column(
    3,
    tags$hr(),
    # fixedPanel(
    wellPanel(
      style = "position:fixed;width:21.5%;",
      tags$h4("Pipeline of TCC-GUI"),
      tags$a("1.  Data input", href = "#Datainput"),
      tags$br(),
      tags$a("2.  Covariance", href = "#Computation"),
      tags$hr(),
      tags$a("3.  Normalization", href = "#MAplot"),
      tags$br(),
      tags$a("4.  Differential Expression", href = "#Volcanoplot"),
      tags$br(),
      tags$a("5.  Enrichment Analysis", href = "#PCAanalysis"),
      tags$br(),
      tags$a("6.  Report", href = "#Heatmap"),
      tags$hr(),
      tags$a("4. More helps", href = "#Morehelps")
    )
  ),
  #column
  column(9,
         tags$hr(),
         tabsetPanel(
           id = "home",
           tabPanel(title = "English", includeMarkdown("document/README_English.md"))
         ))
))
