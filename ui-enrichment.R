#ui-Enrichment analysis

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(column(
  3,
  box(
    title = tagList(icon("cogs"), "Data preparation"),
    id = "DEselect",
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    radioButtons(
      inputId = "filterEN",
      label = "Filter the dataset?",
      choices = c("No", "Yes"),
      selected = "Yes"
    ),
    uiOutput("enrich_filter"),
    do.call(actionBttn, c(
      list(
        inputId = "enrich_sigFilter",
        label = "Run Significance Filter",
        icon = icon("play")),
      actionBttnParams
    )
    )
  ),
  box(
    title = tagList(icon("cogs"), "Parameter Selection"),
    id = "enrich_params",
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    selectInput(
      "orgSelect",
      "Organism",
      choices = c("Mouse" = "org.Mm.eg.db",#mouse
                  "Human Sapien" = "org.Hs.eg.db",#human
                  "Rat" = "org.Rn.eg.db"#rat
                  )
    ),
    selectInput(
      "method",
      "Method",
      choices = c("GO" = "go_enrich",
                  "Kegg" = "kegg_enrich",
                  # "GSEA" = "gsea_enrich",
                  "gseKEGG" = "gseKegg")
    ),
    uiOutput("go_enrich_select"),
    # selectInput(
    #   "keytype",
    #   "Keytype",
    #   choices = c("Pid" = "UNIPROT",
    #               "Gene Name" = "SYMBOL")
    # ),
    numericInput(
      inputId = "enrich_pvalCut",
      label = "p-value cutoff",
      min = 0,
      max = 1,
      value = 0.05,
      step = 0.01,
    ),
    numericInput(
      inputId = "enrich_qvalCut",
      label = "q-value cutoff",
      min = 0,
      max = 1,
      value = 0.05,
      step = 0.01,
    ),
    numericInput(
      inputId = "minGSSize",
      label = "minimal gene set size",
      min = 0,
      max = 10000,
      value = 10,
      step = 1,
    ),
    numericInput(
      inputId = "maxGSSize",
      label = "maximal gene set size",
      min = 0,
      max = 10000,
      value = 800,
      step = 1,
    ),
    do.call(actionBttn, c(
      list(
        inputId = "runEnrichmentAnalysis",
        label = "Run Enrichment Analysis",
        icon = icon("play")),
      actionBttnParams
    )
    )
  )
  
),
column(
  9,
  box(
    title = tagList(icon("table"), "Filtered Result Table"),
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    uiOutput("preEnrichResultTable")
  ),
  box(
    title = tagList(icon("table"), "Enrichment Result Table"),
    width = NULL,
    solidHeader = TRUE,
    status = "primary",
    uiOutput("afterEnrichResultUI")
  ),
  
  tabBox(
    title = "",
    width = NULL,
    tabPanel(title = tagList(icon("sitemap"), "Network"),
             uiOutput("EN_network")),
    tabPanel(
      title = tagList(icon("bar-chart"), "Bubble Plot"),
      uiOutput("EN_dotPlot")
    )
  ))
))