###################
#### Define UIs ###
###################

WGCNA.PPI <- function(){
  fluidRow(
    useShinyjs(),
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
                    } .mypos {
                        margin: 0 auto;
                    }
                    
                    
    ')),
    #imageOutput('poparrow'),
    column(4, class = 'text-center',
           h4(strong(HTML('User Input Specifications'))),
           box(status='primary', background = 'blue', width=40, height=10,
               id='second',
               div(style='display: inline-block;vertical-align:top;',fileInput("input_table", label = strong(HTML("Upload Input File (.xlsx, .csv, .txt)")))),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'FileInfo', label=strong('i'),style='pill',color='default')),
               
               div(style='display: inline-block;vertical-align:top;',selectInput("transformation", h5(strong("Execute Log-2 Transformation of Data?")),
                           choices=list('Yes' = '1',
                                        'No' = '2'),
                           selected = "2")),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TransformInfo', label=strong('i'),style='pill',color='default')),
               
               div(style='display: inline-block;vertical-align:top;',textInput('UserFolder', label=h5(strong("Specify the Name of the Output Folder")), 
                         value='example_output_folder')),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'UserFolderInfo', label=strong('i'),style='pill',color='default'))),
           imageOutput('bufferimagelol'),
           actionButton('ViewExampleFile2', 'View Example File',
                        style='position: absolute; left:10px;background-color:orange')
    ),
    column(4, class = 'text-center',
           uiOutput('ErrorDisplay'),
           h4(strong(HTML('Co-Expression Clustering Parameters'))),
           div(id='divOnTop',uiOutput('exampleFile')),
           div(id='divOnBottom',
             box(status='primary', background = 'blue', width=40, height=10,
                 id='second',
                 div(style='display: inline-block;vertical-align:top;', numericInput("min_cluster_size", label=h5(strong('Minimum Cluster Size')),
                             min=1, max=400, value=10, step=1)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinClusterInfo', label=strong('i'),style='pill',color='default')),
                 div(style='display: inline-block;vertical-align:top;',numericInput("min_kme", label=h5(strong('Minimum kME')),
                             min=0, max=1, value=.7, step=.01)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinKMEInfo', label=strong('i'),style='pill',color='default')),
                 div(style='display: inline-block;vertical-align:top;',numericInput("min_cluster_dist", label=h5(strong('Minimum Cluster Distance')),
                             min=0, max=1, value=.1, step=.01)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinDistInfo', label=strong('i'),style='pill',color='default'))
             ))),
    column(4, class = 'text-center',
           h4(strong(HTML('PPI Analysis Parameters'))),
           box(status='primary', background = 'blue', width=40, height=10,
               id='second',
               div(style='display: inline-block;vertical-align:top;',numericInput("TOM_triggered_module_size", label=h5(strong('Specify How Large a Module Must Be to Trigger TOM')),
                           min=1, max=200, value=100, step=1)),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TomSizeInfo', label=strong('i'),style='pill',color='default')),
               div(style='display:inline-block;vertical-align:top;',actionButton('do', label=h5(strong('Submit WGCNA+PPI Analysis')))),
               div(style='display:inline-block;vertical-align:top;',actionButton('Advanced', 'Advanced Parameters'))),
           imageOutput('ShiftDown'),
           actionButton('ChangeMode', 'Change Analysis Mode',
                        style='display:inline-block;vertical-align:top;position: absolute; right:10px;background-color:orange')
    )
    
    
  )
}

Independent_WGCNA <- function(){
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
      column(4, class = 'text-center',
             h4(strong(HTML('User Input Specifications'))),
             box(status='primary', background = 'blue', width=40, height=10,
                 id='second',
                 div(style='display:inline-block;vertical-align:top;',fileInput("input_table", label = strong(HTML("Upload Input File (.xlsx, .csv, .txt)")))),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'FileInfo', label=strong('i'),style='pill',color='default')),
                 
                 div(style='display:inline-block;vertical-align:top;',selectInput("transformation", h5(strong("Execute Log-2 Transformation of Data?")),
                             choices=list('Yes' = '1',
                                          'No' = '2'),
                             selected = "2")),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TransformInfo', label=strong('i'),style='pill',color='default')),

                 div(style='display:inline-block;vertical-align:top;',textInput('UserFolder', label=h5(strong("Specify the Name of the Output Folder")), 
                           value='example_output_folder')),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'UserFolderInfo', label=strong('i'),style='pill',color='default')),
             
                 ),
             imageOutput('bufferimagelol'),
             actionButton('ViewExampleFile2', 'View Example File',
                          style='position: absolute; left:10px;background-color:orange')
      ),
      column(4,class = 'text-center',
             uiOutput('ErrorDisplay'),
             h4(strong(HTML('Co-Expression Clustering Parameters'))),
             box(status='primary', background = 'blue', width=40, height=10,
                 id='second',
                 div(style='display:inline-block;vertical-align:top;',numericInput("min_cluster_size", label=h5(strong('Minimum Cluster Size')),
                             min=1, max=400, value=10, step=1)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinClusterInfo', label=strong('i'),style='pill',color='default')),
                 
                 div(style='display:inline-block;vertical-align:top;',numericInput("min_kme", label=h5(strong('Minimum kME')),
                             min=0, max=1, value=.7, step=.01)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinKMEInfo', label=strong('i'),style='pill',color='default')),
                 
                 div(style='display:inline-block;vertical-align:top;',numericInput("min_cluster_dist", label=h5(strong('Minimum Cluster Distance')),
                             min=0, max=1, value=.1, step=.01)),
                 div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinDistInfo', label=strong('i'),style='pill',color='default')),
                 div(style='display:inline-block;vertical-align:top;',actionButton('do', label=h5('Submit Independent WGCNA Analysis'))),
                 div(style='display:inline-block;vertical-align:top;',actionButton('Advanced', 'Advanced Parameters'))
             )
      ),
      column(3,
             imageOutput('ShiftDown'),
             actionButton('ChangeMode', 'Change Analysis Mode',
                          style='position: absolute; right:10px;background-color:orange'))
  )
}

Independent_PPI <- function(){
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
    column(4,
           imageOutput('bufferimagelol'),
           actionButton('ViewExampleFile3', 'View Example File',
                        style='position: absolute; left:10px;background-color:orange')),
    column(3, 
    ),
    column(5,class = 'text-center',
           h4(strong(HTML('PPI Analysis Parameters'))),
           box(status='primary', background = 'blue', width=40, height=10,
               id='second',
               div(style='display:inline-block;vertical-align:top;',fileInput("input_table", label = strong(HTML("Upload Input File")))),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'FileInfo', label=strong('i'),style='pill',color='default')),
               
               div(style='display:inline-block;vertical-align:top;',numericInput("TOM_triggered_module_size", label=h5(strong('Specify How Large a Module Must Be to Trigger TOM')),
                                                                                min=1, max=200, value=100, step=1)),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TomSizeInfo', label=strong('i'),style='pill',color='default')),
               div(style='display: inline-block;vertical-align:top;',textInput('UserFolder', label=h5(strong("Specify the Name of the Output Folder")), 
                                                                               value='example_output_folder')),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'UserFolderInfo', label=strong('i'),style='pill',color='default')),
               
               div(style='display:inline-block;vertical-align:top;',actionButton('do', label=h5(strong('Submit Independent PPI Analysis')))),
               div(style='display:inline-block;vertical-align:top;',actionButton('Advanced', 'Advanced Parameters'))
               ),
           imageOutput('ShiftDown'),
           actionButton('ChangeMode', 'Change Analysis Mode',
                        style='position: absolute; right:10px;background-color:orange'))
  )
}

Condensed_UI <- function(){
  fluidRow(
    useShinyjs(),
    column(4, #class = 'text-center',
           uiOutput('ErrorDisplay'),
           h4(strong(HTML('Co-Expression Clustering Parameters'))),
           div(id='divOnTop',uiOutput('exampleFile')),
           div(id='divOnBottom',
               box(status='primary', width=NULL,
                   id='second',
                   div(style='display: inline-block;vertical-align:top;', numericInput("min_cluster_size", label=h5(strong('Minimum Cluster Size')),
                                                                                       min=1, max=400, value=10, step=1)),
                   div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinClusterInfo', label=strong('i'),style='pill',color='default')),
                   div(style='display: inline-block;vertical-align:top;',numericInput("min_kme", label=h5(strong('Minimum kME')),
                                                                                      min=0, max=1, value=.7, step=.01)),
                   div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinKMEInfo', label=strong('i'),style='pill',color='default')),
                   div(style='display: inline-block;vertical-align:top;',numericInput("min_cluster_dist", label=h5(strong('Minimum Cluster Distance')),
                                                                                      min=0, max=1, value=.1, step=.01)),
                   div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'MinDistInfo', label=strong('i'),style='pill',color='default')),
               ),
               #imageOutput('buffeerimageinfinite'),
               div(style='display:inline-block;right:10px;',actionButton('Advanced', 'Advanced Parameters')),
               )),
    column(5,#class = 'text-center',
           h4(strong(HTML('PPI Analysis Parameters'))),
           box(status='primary', width=NULL,
               id='second',
               #div(style='display: inline-block;vertical-align:top;',fileInput("ppi_input_table", label = h5(strong("If executing PPI analysis independently, upload co_exp_clusters_3columns.xlsx file generated by WGCNA or a user-defined file with the same format (.xlsx, .csv, .txt)")))),
               #div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'PPI_FileInfo', label=strong('i'),style='pill',color='default')),
               div(style='display:inline-block;vertical-align:top;',numericInput("TOM_triggered_module_size", label=h5(strong('Maximum Protein Module Size')),
                                                                                 min=5, max=200, value=40, step=1)),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TomSizeInfo', label=strong('i'),style='pill',color='default')),
               div(style='display:inline-block;vertical-align:top;',numericInput("minimum_module_size", label=h5(strong('Minimum Protein Module Size')),
                                                                                 min=1, max=200, value=2, step=1)),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'TomMinInfo', label=strong('i'),style='pill',color='default')),
               div(style='display: inline-block;vertical-align:top;',selectInput("PPIDatabases", h5(strong("Select a PPI Database")),
                                                                                 choices=list('2021 Human Composite BioPlex3+InWeb2016+STRINGv11' = 'BioPlex3_InWeb2016_STRINGv11_human_v1.0.sif',
                                                                                              '2016 Human Composite BioPlex+String400+Inweb150' = 'BioPlex_String400_Inweb150_v1.0.1.sif',
                                                                                              'Human STRINGv11' = 'STRING_v11_human_v1.0.sif',
                                                                                              'Human InWeb2016' = 'InWeb_core2016_v1.0.sif',
                                                                                              'Human BioPlex3' = 'BioPlex3_combined_v1.0.sif',
                                                                                              'Mouse STRINGv11' = 'STRING_v11_mouse_v1.0.sif',
                                                                                              'User-Provided Database' = 'User-Provided Database'),
                                                                                 selected = "BioPlex_String400_Inweb150_v1.0.1.sif")),
               div(style='display: inline-block;vertical-align:top;',actionBttn(inputId = 'DatabaseInfo', label=strong('i'),style='pill',color='default'))
           ),
           #imageOutput('ShiftDown'),
           div(style='display:inline-block;vertical-align:top;position: absolute; right:10px;',
               actionBttn('total_search', 
                          label=h4(strong('Submit JUMPn Analysis')),
                          style='gradient',
                          color='danger'))
    )
  )
}
