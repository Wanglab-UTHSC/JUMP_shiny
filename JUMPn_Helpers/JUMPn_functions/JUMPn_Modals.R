######################################
#### Define Modals for Dynamic UI ####
######################################

inputModal <- function(failed=FALSE){
  modalDialog(class = 'text-center',
    #div(style='display: inline-block;vertical-align:top;',fileInput("meta_file", label = strong(HTML("If your data has more than twenty samples, we recommend uploading a sample-group meta file.")))),
    div(style='display: inline-block;vertical-align:top;',fileInput("meta_file", label = strong(HTML("Upload a meta file")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MetaFileInfo', label=strong('i'),style='pill',color='default')),
    #div(style='display: inline-block;vertical-align:top;',fileInput("ppi_input_file", label = strong(HTML("If running PPI only analysis mode, it is necessary to upload three column PPI input file. See Help  page for formatting.")))),
    div(style='display: inline-block;vertical-align:top;',fileInput("ppi_input_file", label = strong(HTML("Upload co-expression cluster result for 'PPI Only' analysis")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'ppi_input_file_info', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',fileInput("background_file", label = strong(HTML("Upload a background file for pathway enrichment analysis")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'BackgroundFileInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',fileInput("ppi_file", label = strong(HTML("Upload a PPI database")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'PPIFileInfo', label=strong('i'),style='pill',color='default')),
    #div(style='display: inline-block;vertical-align:top;',fileInput("annotation_file", label = strong(HTML("Upload an annotation file for pathway/ontological enrichment")))),
    div(style='display: inline-block;vertical-align:top;',fileInput("annotation_file", label = strong(HTML("Upload an ontology library file for pathway enrichment analysis")))),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'AnnotationFileInfo', label=strong('i'),style='pill',color='default')),
    #div(style='display: inline-block;vertical-align:top;',awesomeRadio('background', label=strong(HTML('Select the background to be used for pathway/ontological enrichment',)),
    div(style='display: inline-block;vertical-align:top;',awesomeRadio('background', label=strong(HTML('Select the background to be used for pathway enrichment analysis',)),
                choices=list('Human/Mouse Background'='In-House Human/Mouse Background (Default)',
                             'User-Supplied Background'='User-Supplied Background',
                             'Merge In-House Background With User-Supplied Background'='Merge In-House Background With User-Supplied Background'),
                selected = c('In-House Human/Mouse Background (Default)'),
                status = 'primary')),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'BackgroundSelectInfo', label=strong('i'),style='pill',color='default')),
    #div(style='display: inline-block;vertical-align:top;',awesomeRadio('PathwayEnrichment', label=strong(HTML('Select databases used for pathway/ontological enrichment')),
    div(style='display: inline-block;vertical-align:top;',awesomeRadio('PathwayEnrichment', label=strong(HTML('Select databases used for pathway enrichment analysis')),
                choices=list('Composite KEGG/GO/Reactome/HALLMARK Database'='Composite KEGG/GO/Reactome/HALLMARK Database',
                             'User-Supplied Database in .xlsx Format'='User'),
                selected=c('Composite KEGG/GO/Reactome/HALLMARK Database'),
                status='primary')),
    div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'AnnotationSelectInfo', label=strong('i'),style='pill',color='default')),
    div(style='display: inline-block;vertical-align:top;',
        awesomeRadio(
      inputId = "ModeofAnalysis",
      label = "Select Mode of Analysis", 
      choices = c("WGCNA+PPI","WGCNA Only" = "Independent WGCNA","PPI Only" = "Independent PPI"),
      selected = "WGCNA+PPI",
      status = "primary"))
  )
}

inputModal3 <- function(failed=FALSE){
  modalDialog(
    h4(strong('DISCLAIMER:',
    'Thank you for choosing JUMPn for your network and pathway analysis.',
    'JUMPn combines Co-Expression clustering with protein-protein interaction data',
    'to readily derive biological insights from complex data. JUMPn users have the choice of',
    'Executing integrated analysis (WGCNA+PPI Analysis; Default) or independent',
    'analysis (Either WGCNA or PPI Analysis separately). Before beginning, please choose between',
    'the available options for analysis.')),
    radioButtons('ModeofAnalysis', label = h5(strong('Select a Mode of Analysis for JUMPn to Execute')),
                 choices=list('WGCNA+PPI'='WGCNA+PPI',
                              'Independent WGCNA'='Independent WGCNA',
                              'Independent PPI'='Independent PPI')),
    footer = actionButton('MODA_Begin', 'Begin')
  )
}
testSweetAlert <- function(session){
  sendSweetAlert(
    session=session,
    title='Information',
    type='info',
    btn_labels = NA,
    text = tags$div(class='text-center','Thank you for choosing JUMPn for your network and pathway analysis.JUMPn combines Co-Expression clustering with protein-protein interaction data
            to readily derive biological insights from complex data. JUMPn users have the choice of executing integrated analysis (WGCNA+PPI Analysis; Default) or independent
            analysis (Either WGCNA or PPI Analysis separately). Before beginning, please choose between the available options for analysis.',
            awesomeRadio('ModeOfAnalysis', label = h5(strong('Select a Mode of Analysis for JUMPn to Execute')),
                          choices=list('WGCNA+PPI'='WGCNA+PPI',
                                        'Independent WGCNA'='Independent WGCNA',
                                        'Independent PPI'='Independent PPI'),
                         width='100%'),
            div(style='display:inline-block;vertical-align:top;',actionButton('MODA_Begin', 'Begin'))
            
    )
    
  )
}

messageModal <- function(type,session){
  sendSweetAlert(
    session=session,
    type='info',
    title='Information',
    btn_labels = NA,
    text = tags$div(
      h4(strong(sprintf('Great, you can begin the %s mode of analysis. Check the proper file format with the "View Example File" button below', type))),
      actionButton('ViewExampleFile', 'View Example File'),
      actionButton('Continue', 'Continue')
      
    )
    
  )
  
}

imageModal <- function(output){
  modalDialog(
    #h4(strong("For WGCNA+PPI analysis or WGCNA Only analysis, the Gene List (Gene Symbols) must be in Column 2. Expression data must be in Column 4 onward; you can use as many expression columns as desired. For analyses of more than 20 samples, we recommend providing a group meta file as described in the 'Advanced Parameters' panel.")),
    h4(strong("Column 1: unique feature ID; column 2: official gene symbol that match the ontology library and PPI networks; column 3: can be blank or flexible to users; from column 4 onwards: protein abundance of individual samples (See an EXAMPLE of file format in the attached image below). 
	The number of samples is not limited by the program but may depend on computer memory (4Gb good for sample size < 30; 16Gb for n = 200 samples). 
	For analyses of large sample size dataset (e.g. >30 samples), we recommend providing a group meta file as described in the 'Advanced Parameters' panel.")),
    renderImage({
      list(src='www/images/WGCNAInput.png',
            contentType = 'image/png',
            width = 750,
            height = 600)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}



imageModal2 <- function(output){
  modalDialog(
    h4(strong("PPI Analysis requires a different input than is needed for WGCNA. Column 1 must have the Cluster/Group to which the gene belongs. ",
              "Column 2 must have the accession ID (or some other misscelaneous categorical information), and Column 3 must contain the gene list (Gene Symbols),",
              'If you are conducting independent PPI analysis after running WGCNA, then this 3-column file is produced automatically, which you  can',
              'use directly as input for independent PPI after downloading the WGCNA data. Otherwise, you can make your own custom file which adheres to the described format.',
              'Additional columns may be provided beyond the standard 3-column format,  but these columns will not be evaluated.')),
    renderImage({
      list(src='www/images/PPIInput.png',
           contentType = 'image/png',
           width = 550,
           height = 425)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}
imageModal3 <- function(output){
  modalDialog(
    h4(strong("A user-definede PPI database must be a .sif file with a format as in the image below. Column 1 contains interaction component A
              Column 2 contains any misscelaneous information, and Column 3 contains interaction component B.")),
    renderImage({
      list(src='www/images/example_ppi_input.png',
           contentType = 'image/png',
           width = 450,
           height = 600)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}
imageModal4 <- function(output){
  modalDialog(
    h4(strong("A user-definede pathway annotation database must have the following format. Column 1 contains the pathway name or description,
              Column 2 contains the genes in the pathway separated by '/'")),
    renderImage({
      list(src='www/images/example_pathway_enrichment.png',
           contentType = 'image/png',
           width = 650,
           height = 400)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}
imageModal5 <- function(output){
  modalDialog(
    h4(strong("A user-definede background file must be a .txt file with simply a list of genes. It should have a format as appears in the image below.")),
    renderImage({
      list(src='www/images/example_background.png',
           contentType = 'image/png',
           width = 650,
           height = 400)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}

imageModal6 <- function(output){
  modalDialog(
    h4(strong("A meta-group file must have the identical sample names from the input expression matrix in Column 1. Column 2 onwards can be any group identifying information for each sample.")),
    renderImage({
      list(src='www/images/example_meta_file.png',
           contentType = 'image/png',
           width = 650,
           height = 400)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}

imageModal7 <- function(output){
  modalDialog(
    h4(strong("These are the file types allowed by JUMPn.")),
    renderImage({
      list(src='www/images/acceptable_file_formats.png',
           contentType = 'image/png',
           width = 650,
           height = 300)
    },deleteFile=FALSE),
    easyClose = TRUE
  )
  
}
programFinishedModal <- function(session, location){
  sendSweetAlert(
    session = session,
    title='Success !!', 
    type='success',
    btn_labels = NA,
    text = tags$div(sprintf('Results have been  generated. The data has been stored at: \n %s ', location),
                    actionButton('TakeMeThere', 'Continue to Results')
    )
  )
}
NewPPISelected <- function(session){
  sendSweetAlert(
    session = session,
    title='Warning!', 
    type='warning',
    btn_labels = NA,
    text = tags$div("PPI Analysis requires a different input than is needed for WGCNA. Column 1 must have the Cluster/Group to which the gene belongs. ",
              "Column 2 must have the accession ID (or some other misscelaneous categorical information), and Column 3 must contain the gene list (Gene Symbols),",
              'If you are conducting independent PPI analysis after running WGCNA, then this 3-column file is produced automatically (co_exp_clusters_3columns.xlsx), which you  can',
              'use directly as input for independent PPI after downloading the WGCNA data. Otherwise, you can make your own custom file which adheres to the described format.',
              'Visit the "Help" page for an example input PPI file.'
    ),
    
    easyClose=TRUE
  )
}

BCellProtModal <- function(session,wd){
  sendSweetAlert(
    session = session,
    title='Success !!', 
    type='success',
    btn_labels = NA,
    text = tags$div(sprintf('The demo B-Cell proteomic file %s has been uploaded. You may now adjust the search parameters and submit the analysis.',wd)
    ),
    
    easyClose=TRUE
  )
}

PPISelected <- function(output){
  modalDialog(
    h4(strong('You selected the independent PPI mode of analysis. Make sure you upload a three column input file with Column 1 having the group information for each protein,
                     Column 2 has miscellaneous information, and Column 3 has the gene/protein names. A file of this format is automaticlaly generated by WGCNA (co_exp_clusters_3columns.xlsx)')),
    renderImage({
      list(src='www/images/PPIInput.png',
           contentType = 'image/png',
           width = 450,
           height = 600)
    }, deleteFile=FALSE),
    easyClose=TRUE
  )
}

ArrivedToResultsModal <- function(whereami,session){
  if (whereami=='WGCNA+PPI'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your WGCNA results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
          Once you have finished investigating your results, you can navigate to the PPI results via the Page 2 tab above.',
        actionButton('ViewResults', 'View Results')
      ),
      easyClose=TRUE
    )
  } else if (whereami=='Independent WGCNA'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your WGCNA results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
          Once you have finished investigating your results, you can submit another search  by navigating back to the  Commence Analysis Page.',
        actionButton('ViewResults', 'View Results'))
      )
  } else if (whereami=='Independent PPI'){
    sendSweetAlert(
      session=session,
      type='info',
      title='Information',
      btn_labels = NA,
      text = tags$div(
        'Here are your PPI results. Below you will find co-expression trends (in the form of boxplot and linear figures), as well as pathway analysis (if selected). 
        Once you have finished investigating your results, you can submit another search  by navigating back to the  Commence Analysis Page.'
      )
    )
  }
}

QuickPPIModal <- function(){
  fluidRow(
    tags$style(HTML('#btn1 {
                        background-color:red
                    } #btn2 {
                        background-color:red
                    } #Advanced {
                        background-color:lightblue
                    } #PPISearch {
                        background-color:lightblue
                    } #WGCNASearch {
                        background-color:lightblue
                    } #do {
                        background-color:red
                    
                    
    ')),
    column(4),
    column(4, 
           h4(strong(HTML('PPI Analysis Parameters'))),
           box(status='primary', background = 'blue', width=40, height=10,
               id='second',
               fileInput("input_table", label = strong(HTML("Upload Input File"))),
               selectInput("skip_PPI_network_analysis", h5(strong("Should Co-Expression Clusters Be Modularized by Protein Interactions?")),
                           choices=list('Yes'='Yes',
                                        'No'='No')),
               sliderInput("TOM_triggered_module_size", label=h5(strong('Specify How Large a Module Must Be to Trigger TOM')),
                           min=1, max=200, value=100, step=1, ticks=TRUE),
               div(style='display:inline-block;vertical-align:top;',actionButton('Advanced', 'Advanced Parameters')),
               div(style='display:inline-block;vertical-align:top;',actionButton('do', label=h4(strong('Submit Independent PPI Analysis')))))
           
    ),
    column(4)
  )
}
ConfirmSearchModal <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    text = tags$div('You are about to submit a search. Is your file of the correct format and are your parameters set. If so click "Confirm Search". Otherwise press "Back" and double check your parameters.',
                    actionButton("ConfirmSearch", "Confirm Search"),
                    actionButton('BackButton', 'Back'))
  )
}

Upload_error <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    text = tags$div('Whoops! It appears you  forgot to upload an input file',
                    actionButton('BackButton2', 'Back'))
  )
}

ConfirmNewSearchModal <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    text = tags$div('Submiting a new search will cause data from the previous run to be  overwritten. Wee highly recommend, if you have not already, to download the data generated with the
                    download handler below. Otherwise, you can confirm your new search with the Continue with New Search button',
                    actionButton("ConfirmNewSearch", "Confirm New Search"),
                    actionButton('CancelNewSearch', 'Cancel'))
  )
}

WidgetInfoModal <- function(session,messageinfo){
  sendSweetAlert(
    session=session,
    title='Information',
    text=messageinfo,
    type='info',
  )
}
ChangeFolder <- function(session){
  sendSweetAlert(
    session=session,
    type = 'warning',
    title='Warning:',
    btn_labels = NA,
    showCloseButton = TRUE,
    text = tags$div('We detected an output folder with the identical name currently specified. If you do not  wish for this data to be overwritten, please change it here, or simply continue',
                    textInput('UserFolder', label=h5(strong("Specify the Name of the Output Folder for Data Export")), 
                              value='example_output_folder'),
                    actionButton('total_search', 'Continue Search'))
  )
}

  
successful_download <- function(session){
  sendSweetAlert(
    session=session,
    title = 'Success',
    text='Congrats, publciation data has been downloaded!',
    type='success'
  )
}


  
