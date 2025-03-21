library(shiny)
library(shinythemes)
library(fresh)

#Acknowledgement: 
#Some code of this program was adapted from a program related under the MIT license:
#TCC-GUI(https://github.com/swsoyee/TCC-GUI)

tagList(
  dashboardPage(
    header = dashboardHeader(
        tags$li(
          class = "dropdown",
          tags$style(".main-header {min-height: 50px}"),
        ),
        title = "JUMP Shiny"
    ),
    dashboardSidebar(
      sidebarMenu(
        id = "sider",
        menuItem("Home",
                 tabName = "welcome",
                 icon = icon("home")
        ),
        menuItem("Tutorial",
                 tabName = "manual",
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
          "Differential Expresssion",
          tabName = "diffExpression",
          icon = icon("filter")
        ),
        menuItem(
          "Enrichment Analysis",
          tabName = "enrichmentMethods",
          icon = icon("circle-nodes")
        )

      )
    ),
    dashboardBody(
      includeCSS("www/style/theme.css"),
      tabItems(
        tabItem(
          tabName = "welcome",
          tabBox(
            title = "",
            width = NULL,
              fluidRow(
                column(
                  includeMarkdown("document/English_Welcome.md"),
                  width = 10,
                  offset = 1
                )
              )
          )
        ),
        tabItem(
          tabName = "manual",
          tabBox(
            title = "",
            width = NULL,
            tabPanel(
              title = "Experiment Design",
              icon = icon("pen-to-square"),
              fluidRow(column(
                includeMarkdown("document/English_Experiment_design.md"),
                width = 10,
                offset = 1
              ))
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
              title = "Differential Expression",
              icon = icon("filter"),
              fluidRow(column(
                includeMarkdown("document/English_Expression.md"),
                width = 10,
                offset = 1
              ))
            ),
            tabPanel(
              title = "Enrichment pathway analysis",
              icon = icon("circle-nodes",lib = "font-awesome"),
              fluidRow(column(
                includeMarkdown("document/English_Enrichment.md"),
                width = 10,
                offset = 1
              ))
            )
          )
        ),
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
        )$value)

      )
    )
  )

)
