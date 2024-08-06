# ui-report.R

fluidPage(
  includeCSS("www/style/theme.css"),
  column(
  3,
  box(
    title = tagList(icon("cogs"), "Report Parameters"),
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    footer = "HTML report (default) is highly recommended.",
    radioGroupButtons(
      inputId = 'format',
      label = 'Document Format',
      choices = c('HTML','docx'),
      justified = TRUE,
      status = "primary"
    ),
    do.call(actionBttn, c(
      list(
        inputId = "generateReport",
        label = "Generate Report",
        icon = icon("play")
      ),
      actionBttnParams
    )),
    
    uiOutput("renderDownloadButton"),
    uiOutput("renderReportButton")
  )
),
column(9,
       box(
         title = "Output Option",
         status = "primary",
         solidHeader = TRUE,
         width = NULL,
         textInput(
           inputId = "projName",
           label = "Project Name",
           value = "Proteomics Analysis",
           placeholder = "Enter your project name."
         ),
         textInput(
           inputId = "speciesN",
           label = "Species",
           value = "Homo Sapiens",
           placeholder = "Enter species of data."
         ),
         textInput(
           inputId = "tissueType",
           label = "Tissue Type",
           value = "Brain",
           placeholder = "Enter tissue type."
         )
       )))