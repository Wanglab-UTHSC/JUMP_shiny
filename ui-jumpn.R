library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyWidgets)
require(visNetwork, quietly = TRUE)
options(shiny.maxRequestSize = 30*1024^2)
R_BUILD_TAR=tar

# Define UI for application that draws a histogram
ui <- fluidPage(
    h1(strong(HTML('Welcome to JUMPn Network Analysis'))),
    #titlePanel("Welcome to JUMPn Network Analysis"),
    fluid = TRUE,
    collapsible=TRUE,
    
    tabsetPanel(
        id='whichtab',

        tabPanel(
            "Parameter settings", 
            value='Commence',
            #includeHTML("JUMPn_Helpers/html_files/network_analysis.html"),
            h2(strong(HTML(''))),
            uiOutput('DynamicUI')),
        
        tabPanel(
            title=uiOutput('WGCNATitlePage'),
            value='WGCNAPage',
            #includeHTML("JUMPn_Helpers/html_files/net_output.html"),
            useShinyjs(),
            inlineCSS(list(.networkplotborder1 = "border: 2px solid black")),
            inlineCSS(list(.networkplotborder2 = "border: 2px solid black")),
            mainPanel(
                uiOutput('ResultsPath'),
                uiOutput('ExpressionDisplayOptions'),
                
                column(8,
                       uiOutput('ExpressionDisplayObject')),
                column(4,
                       uiOutput('ClusterPathwayHeatmap')),
                
                imageOutput('buffeerex'),
                imageOutput('buffeerey'),
                
                uiOutput('ClusterDTSelectionBox'),
                uiOutput('DynamicClusterDTDisplay'),
                fluidRow(
                  column(6,
                         plotOutput('GeneExpression')),
                  column(6,
                         plotOutput('Cluster_Expression'))
                ),
                htmlOutput('ClusterPathwayTableTitle'),
                dataTableOutput('ClusterPathwayTable'),
                
                
            )
        ),
        
        tabPanel(
            title=uiOutput('PPITitlePage'),
            value='PPIPage',
            #includeHTML("JUMPn_Helpers/html_files/path_output.html"),
            mainPanel(
              uiOutput('PPIResultsPath'),
              fluidRow(
                column(8,
                  uiOutput('DynamicClusterSelectionBox')),
                column(4,
                  uiOutput('TheClusterExpressionDisplay'),
                )
              ),
              fluidRow(
                column(8,
                       div(style='display:inline-block;vertical-align:top;',uiOutput('DynamicClusterDisplay'))),
                                                                                         
                column(4,
                       uiOutput('ExpressionFormat'))
              ),

              fluidRow(
                column(8,
                  uiOutput('DynamicModuleSelectionBox')),
                column(4,
                  uiOutput('ModulePathwaySelection'))
              ),
              fluidRow(
                column(8,
                       uiOutput('DynamicModularDisplay')),
                column(4,
                       uiOutput('EMAP'))
              )
            )
        ),
    )
)
