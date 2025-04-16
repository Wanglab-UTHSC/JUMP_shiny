fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(
    column(
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
        radioButtons(
          inputId = "ImputationSelection",
          label = "Missing Value Imputation",
          choices = c(
            "No Imputation" = "NA",
            "Imputation" = "Imputation",
            "Data without NAs" = "CleanData"
          ),
          selected = "NA"
        ),
        conditionalPanel(
          condition = "input.ImputationSelection == 'Imputation'",
          numericInput(
            inputId = "filterValue",
            label = "Enter Filter Value",
            value = 1
          )
        ),
        # Popovers for each radio button choice
        radioTooltip("ImputationSelection", choice = "NA", title = "No imputation on whole data"),
        radioTooltip("ImputationSelection", choice = "Imputation", title = "See documentation for details"),
        radioTooltip("ImputationSelection", choice = "CleanData", title = "No imputation and perform DE only on data without NAs"),
        do.call(
          actionBttn, c(
            list(
              inputId = "DETestType",
              label = "Run Differential Expression Test",
              icon = icon("play")
            ),
            actionBttnParams
          )
        )
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
      tabBox(
        title = "",
        width = NULL,
        # tabPanel(
        #   title = tagList(icon("bar-chart"), "Moving SD"),
        #   uiOutput("DE_distributionUI")
        # ),
        # tabPanel(
        #   title = tagList(icon("bar-chart"), "SD within group"),
        #   uiOutput("DE_OldSDUI")
        # ),
        tabPanel(
          title = tagList(icon("bar-chart"), "Volcano Plot"),
          fluidRow(
            column(
              3,
              textInput("volplttitle", "Graphic Title", value = "Volcano Plot"),
              selectInput(
                "metric2",
                label = "Select the measure of significance",
                choices = c("p-value" = "p-value", "FDR" = "FDR"),
                selected = "p-value"
              ),
              numericInput(
                "pvalcutoff",
                label = "Significance level",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.05
              ),
              numericInput(
                "log2fcUp",
                label = "Upregulated Log2-fold cutoff",
                min = -2,
                max = 2,
                step = 0.1,
                value = 0.1
              ),
              numericInput(
                "log2fcDown",
                label = "Downregulated Log2-fold cutoff",
                min = -2,
                max = 2,
                step = 0.1,
                value = -0.1
              ),
              do.call(
                actionBttn, c(
                  list(
                    inputId = "RunVolcano",
                    label = "Run Volcano Plot",
                    icon = icon("play")
                  ),
                  actionBttnParams
                )
              ),
              box(
                title = tagList(icon("option-vertical",lib = "glyphicon"), "Advanced Parameters"),
                solidHeader = F,
                width = NULL,
                collapsible = T,
                collapsed = T,
                textAreaInput(
                  inputId = "labelProt",
                  label = "Interested protein to label on plot",
                  rows = 5,
                  placeholder = "Input your protein accession number/Uniprot ID.\nOne protein each line"
                ),
                textInput(inputId = "volcxlab", 
                          "x-axis Label", 
                          value = "log<sub>2</sub>"),
                sliderInput(
                  "volcanoPointSize",
                  "Point Size",
                  min = 1,
                  max = 8,
                  value = 3,
                  step = 0.2
                ),
                spectrumInput(
                  inputId = "downColor",
                  label = tagList("Down-regulated in the second group", htmlOutput("downPreview")),
                  choices = list(
                    list(
                      "#0000ff",
                      "green",
                      "black",
                      "grey30",
                      "white"
                    ),
                    as.list(brewer.pal(n = 9, name = "Blues")),
                    as.list(brewer.pal(n = 9, name = "Greens")),
                    as.list(brewer.pal(n = 9, name = "Greys"))
                  ),
                  options = list(`toggle-palette-more-text` = "Show more")
                ),
                spectrumInput(
                  inputId = "upColor",
                  label = tagList("Up-regulated in the second group", htmlOutput("upPreview")),
                  choices = list(
                    list(
                      "red",
                      "yellow",
                      "orange",
                      "white",
                      "#0000ff"
                    ),
                    as.list(brewer.pal(n = 9, name = "Reds")),
                    as.list(brewer.pal(n = 9, name = "Oranges")),
                    as.list(brewer.pal(n = 9, name = "Greys"))
                  ),
                  options = list(`toggle-palette-more-text` = "Show more")
                )
              )
            ),
            column(
              9,
              downloadButton("downloadVolcanoTable", "Download DE Table"),
              plotlyOutput("volcanoPlot") %>% withSpinner(),
              tags$hr(),
              plotlyOutput("protBarPlotInVolcano") %>% withSpinner()
            )
          )
        ),
        tabPanel(
          title = tagList(icon("sitemap"), "Heatmap"),
          tagList(fluidRow(
            column(
              3,
              numericInput("heatmap_log2fc", 
                           label = "Log2-fold cutoff", 
                           min = 0, 
                           max = 10, 
                           step = 0.1,
                           value = 0),
              textInput(
                inputId = "heatmapTitle",
                label = "Heatmap Title",
                value = "Heatmap of top proteins"
              ),
              checkboxInput(
                inputId = "showrowlabel",
                label = "Show Row Label?",
                value = F),
              checkboxInput(
                inputId = "rowlabelUpper",
                label = "Upper-case row labels?",
                value = F),
              checkboxInput(
                inputId = "showdendrow",
                label = "Show Row Dendrogram?",
                value = T),
              checkboxInput(
                inputId = "hm_NA",
                label = "Remove protein/peptides with NA values?",
                value = F
              ),
              do.call(actionBttn, c(
                list(
                  inputId = "heatmapRun",
                  label = "Run Heatmap",
                  icon = icon("play")
                ),
                actionBttnParams
              ))
            ),
            column(
              9,
              downloadButton("downloadHeatmap", "Download Heatmap Table"),
              plotlyOutput("heatmap")%>% withSpinner()
            )
          ))
        )
      )
    )
  )
)

