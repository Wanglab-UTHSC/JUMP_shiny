#ui-Enrichment analysis

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(column(
  3,
  box(
    title = tagList(icon("sliders-h"), "Data source"),
    width = NULL, solidHeader = TRUE, status = "primary",
    selectInput(
      "enrichGeneSource",
      "Genes to use:",
      choices = c("Differentially expressed (from DE result above)" = "de",
                  "Custom list (paste)"= "custom"),
      selected = "de"
    ),
    conditionalPanel(
      condition = "input.enrichGeneSource == 'de'",
      tagList(radioButtons(
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
          icon = icon("play")
        ),
        actionBttnParams
      )))
    ),
    conditionalPanel(
      condition = "input.enrichGeneSource == 'custom'",
      tagList(
        textAreaInput(
          "customGenes",
          "Paste gene symbols/Uniprot ID/Ensembl/RefSeq (one per line, or comma/space/semicolon separated):",
          placeholder = "CNR1, BRCA1, APOE MAPT; mTOR",
          rows = 8
        ),
        do.call(actionBttn, c(
          list(
            inputId = "useCustomGenes",
            label   = "Use this list for enrichment",
            icon    = icon("check")
          ),
          actionBttnParams
        ))
      )
    )
  ),
  box(
    title = tagList(icon("database"), "Custom Pathway Database (optional)"),
    width = NULL, solidHeader = TRUE, status = "info",
    fileInput(
      inputId = "customPathwayFile",
      label = "Upload custom pathway database:",
      accept = c(".csv", ".txt", ".tsv"),
      placeholder = "Accept plain txt/csv file...",
      buttonLabel = "Browse...",
      multiple = FALSE
    ),
    
    radioButtons(
      "pathwayFileFormat",
      "File format:",
      choices = c("Standard(Pathway, Description, Gene)" = "std",
                  "Two-column(Pathway, Gene)" = "two_column",
                  "Multi-column(Pathway in col1, genes in other cols)" = "multi_column"),
      selected = "std"
    ),
    
    tags$div(
      style = "font-size: 12px; color: #gray; margin-top: -5px; margin-bottom: 10px;",
      conditionalPanel(
        condition = "input.pathwayFileFormat == 'std'",
        "Standard format: Pathway,Description,Gene (one gene per row)"
      ),
      conditionalPanel(
        condition = "input.pathwayFileFormat == 'two_column'",
        "Two columns: Pathway,Gene (one gene per row)"
      ),
      conditionalPanel(
        condition = "input.pathwayFileFormat == 'multi_column'",
        "Multiple columns: Pathway,Gene1,Gene2,Gene3 ..."
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
                  "Custom (uploaded)" = "custom_enrich")
    ),
    uiOutput("go_enrich_select"),
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