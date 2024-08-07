library(shiny)
library(shiny)
library(shinythemes)
library(fresh)

stjudelogo <- tags$a(href='https://www.stjude.org',
  tags$img(
  src = "images/icons/SJ_Full_H_W.png",
  style = 'height: 50px; width: 120px; left: 50%;position: absolute; transform: translateX(-50%);'
),target = "_blank")



header <-  htmltools::tagQuery(dashboardHeader(
  tags$li(
    class = "dropdown",
    tags$style(".main-header {min-height: 50px}"),
  ),
  title = "JUMP Shiny"))
header <- header$
  addAttrs(style = 'position:relative;')$ # add some styles to the header 
  find(".navbar.navbar-static-top")$ # find the header right side
  append(stjudelogo)$ # 
  allTags()

tagList(
  dashboardPage(
    header = header,
    dashboardSidebar(
      sidebarMenu(
        id = "sider",
        menuItem("Documentation",
                 tabName = "introduction",
                 icon = icon("book")
        ),
        menuItem(
          "Experiment Design",
          tabName = "BlockRand",
          icon = icon("pen-to-square")
        ),
        menuItem(
          "Exploratory Analysis",
          tabName = "dateImport",
          icon = icon("file-import")
        ),
        menuItem(
          "Batch Normalization",
          tabName = "normalizationMethods",
          icon = icon("align-justify")
        ),
        menuItem(
          "Covariate Analysis",
          tabName = "covTab",
          icon = icon("chart-line")
        ),
        menuItem(
          "Differential Expresssion",
          tabName = "diffExpression",
          icon = icon("filter")
        ),
        menuItem(
          "Enrichment Method",
          tabName = "enrichmentMethods",
          icon = icon("circle-nodes")
        ),
        # menuItem(
        #   "Network Analysis",
        #   tabName = "NetworkAnalysis",
        #   icon = icon("network-wired")
        # ),
        menuItem(
          "Report Export",
          tabName = "reportTab",
          icon = icon("download")
        )
      )
    ),
    dashboardBody(
      includeCSS("www/style/theme.css"),
      tabItems(
        tabItem(
          tabName = "introduction",
          tabBox(
            title = "",
            width = NULL,
            tabPanel(
              title = "Welcome to JUMP Shiny",
              icon = icon("info"),
              fluidRow(
                column(
                  includeMarkdown("document/English_Welcome.md"),
                  width = 10,
                  offset = 1
                )
              )
            ),
            tabPanel(
              title = "Exploratory Analysis",
              icon = icon("file-import"),
              fluidRow(column(
                includeMarkdown("document/English_Data_input.md"),
                width = 10,
                offset = 1
              ))
            ),
            tabPanel(
              title = "Batch Normalization",
              icon = icon("align-justify"),
              fluidRow(column(
                includeMarkdown("document/English_Normalization.md"),
                width = 10,
                offset = 1
              ))
            ),
            tabPanel(
              title = "Covariate Analysis",
              icon = icon("chart-line"),
              fluidRow(column(
                includeMarkdown("document/English_Covariance.md"),
                width = 10,
                offset = 1
              ))
            ),
            # tabPanel(
            #   title = "eSEM Analysis",
            #   icon = icon("magnifying-glass-chart"),
            #   fluidRow(column(
            #     includeMarkdown("document/English_eSEM.md"),
            #     width = 10,
            #     offset = 1
            #   ))
            # ),
            tabPanel(
              title = "Differential Expression",
              icon = icon("filter"),
              fluidRow(column(
                includeMarkdown("document/English_Expression.md"),
                width = 10,
                offset = 1
              ))
            ),
            tabPanel(
              title = "Enrichment Method",
              icon = icon("circle-nodes"),
              fluidRow(column(
                includeMarkdown("document/English_Enrichment.md"),
                width = 10,
                offset = 1
              ))
            )
            # tabPanel(
            #   title = "Report Export",
            #   icon = icon("download"),
            #   fluidRow(
            #     column(
            #       includeMarkdown("document/English_Report.md"),
            #       width = 10,
            #       offset = 1
            #     )
            #   )
            # )
          )
        ),
        tabItem(tabName = "guidence", source(
          file = "ui-homepage.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "BlockRand", source(
          file = "ui-blockRand.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "dateImport", source(
          file = "ui-data-import.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "normalizationMethods", source(
          file = "ui-normalization.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "diffExpression", source(
          file = "ui-differential-expression.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "enrichmentMethods", source(
          file = "ui-enrichment.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        # tabItem(tabName = "NetworkAnalysis", source(
        #   file = "ui-jumpn.R",
        #   local = TRUE,
        #   encoding = "UTF-8"
        # )$value),
        tabItem(tabName = "covTab", source(
          file = "ui-covariance.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value),
        tabItem(tabName = "reportTab", source(
          file = "ui-report.R",
          local = TRUE,
          encoding = "UTF-8"
        )$value)
      )
    )
  ),
  tags$footer(
    class = "main-footer",
    div(
      class = "footer-links",
      tags$a(href = "https://www.stjude.org/legal/st-jude-privacy-policy-statement.html", "U.S. Privacy Notice"),
      tags$a(href = "https://www.stjude.org/legal.html", "Disclaimer / Registrations / Copyright Statement"),
      tags$a(href = "https://www.stjude.org/legal/linking-policy.html", "Linking Policy"),
      tags$a(href = "https://www.stjude.org/legal/notice-of-privacy-practices.html", "Notice of Privacy Practices (HIPAA)"),
      tags$a(href = "https://www.stjude.org/legal/notice-of-non-discrimination.html", "Notice of Non-Discrimination")
    ),
    div(
      "© Copyright 2024. St. Jude Children's Research Hospital, a not-for-profit, section 501(c)(3)."
    )
  )
)
