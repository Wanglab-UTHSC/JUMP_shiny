
########################################################################
## JUMPn Script for Error Handling With Sophisticated Alerts to User ###
########################################################################

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

Check_file_type <- function(session, input_file, acceptable_format_vector, jump_file, check_vector=NA){
  extension <- file_ext(input_file)
  acceptable_extension <- paste(acceptable_format_vector, sep='', collapse = ',')
  print(acceptable_extension)
  print(extension)
  error_ctr <- 0
  if (is.element(extension,acceptable_format_vector)==FALSE){
    sendSweetAlert(
      session=session,
      type = 'warning',
      title='Warning:',
      btn_labels = NA,
      text = tags$div(sprintf('Whoops! It appears that file is not of one of the acceptable types: %s',acceptable_extension)),
      closeOnClickOutside = TRUE
    )
    error_ctr <-error_ctr + 1
  } else {
    if (extension == 'xlsx'){                                                                            ### If the User File is .xlsx ...
      user_db <- read_excel(input_file)                                                               ### Call read_excel function on the user input file
    } else if (extension == 'csv'){                                                                     ### If the User file is .csv  ...
      user_db <- read.csv(input_file)                                                              ### Also call read_excel function
    } else if (extension == 'txt'){                                                                     ### If the User file is ,txt
      print('1')
      user_db <- read_tsv(input_file)                                                               ### Call base r read.table function
      print('2')
    } else if (extension == 'sif'){                                                                     ### If the User file is ,txt
      user_db <- read.table(input_file)                                                               ### Call base r read.table function
    }
    summary(user_db)
    if (jump_file == 'wgcna'){
      ##############################################################
      ###### Checking to see if user should upload  meta file ###### 
      ##############################################################
      if (ncol(user_db) >24){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('It appears your input file has more than 20 samples. We recommend when doing large batch analysis to include a meta file containing group info. This is optional, but if you have a meta file available, go to the Advanced Parameters panel to upload the file or go to the Help page to view the file parameters')),
          closeOnClickOutside = TRUE
        )
      }
      ###############################################################
      ###### Checking to see if columns 4 onward are numerical ###### 
      ###############################################################
      if (ncol(user_db) <  4){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('It appears your file has less than four columns. Are you accidentally uploading a three column PPI file to the WGCNA file widget? If so, please upload the three column file in the Advanced Parameters panel. Otherwise, consult the Help page for more information on file formating.')),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      } else {
      expression_data <- as.data.frame(user_db[,4:ncol(user_db)])
      if (ncol(expression_data)<4){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('JUMPn requires a minimum of four expression columns in order to execute WGCNA. Please add more expression replicates.')),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      }
      #for (col in 1:ncol(expression_data)){
       # print(is.numeric(expression_data[1,col]))
        #if (is.numeric(expression_data[1,col])==FALSE){
         # print(expression_data[2,col])
          #sendSweetAlert(
           # session=session,
            #type = 'warning',
            #title='Warning:',
            #btn_labels = NA,
            #text = tags$div(sprintf('Columns 4 onward should be numerical data columns, but it appears these columns contain text in your file. Please consult the help page to correct the file format')),
            #closeOnClickOutside = TRUE
          #)
          #error_ctr <-error_ctr + 1
        #}
      #}
      print(length(intersect(toupper(unlist(user_db[,2])),check_vector))/nrow(expression_data))
      print(length(intersect(toupper(unlist(user_db[,2])),check_vector)))
      print(nrow(expression_data))
      print(ncol(expression_data))
      print(expression_data)
      if (length(intersect(toupper(unlist(user_db[,2])),check_vector))/nrow(expression_data)<.05) {
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('Woah. Your input genes matched with less than 5 percent of our in-house background. This may yield poor hits during pathway analysis? Is your data coming from a species other than human or mouse? If so, we recommend uploading your own backgroun and pathway database.')),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      }
      }
    } else if (jump_file == 'meta_file'){
      if (length(check_vector)>=8){
      sample_names <- check_vector[4:length(check_vector)]
      meta_names <- as.vector(unlist(user_db[,1]))
      #print(user_db)
      if (ncol(user_db) < 2){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('Whoops. It appears your meta file only has one column. Did you only provide sample names? JUMPn requires that sample names be associated with at least one group when conducting meta analysis.')),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      } else {
      if (length(sample_names) != length(meta_names)){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('Whoops. It appears your input file contains %d samples while your meta file contains %d samples. These files must have identical sample numbers with identical sample names. This must be fixed in order to continue', length(sample_names),length(meta_names))),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      } else {
        if (sum(sample_names==meta_names)!=nrow(user_db)) {
          sendSweetAlert(
            session=session,
            type = 'warning',
            title='Warning:',
            btn_labels = NA,
            text = tags$div(sprintf('Whoops. It appears the sample names in your meta file do not match the names in your input matrix. You will need to fix this to continue')),
            closeOnClickOutside = TRUE
          )
          error_ctr <-error_ctr + 1
        }
      }
      }
      } else {
      sendSweetAlert(
        session=session,
        type = 'warning',
        title='Warning:',
        btn_labels = NA,
        text = tags$div(sprintf('Whoops. It appears your input quantification file still does not have the correct format. Either your file does not have the minimum 8 total columns, does not have the minimum four expression replicates, or the expression data is not entirely numerical.')),
        closeOnClickOutside = TRUE
      )
      error_ctr <-error_ctr + 1
    }
  } else if (jump_file == 'ppi_input_file'){
    if (ncol(user_db)!=3){
      sendSweetAlert(
        session=session,
        type = 'warning',
        title='Warning:',
        btn_labels = NA,
        text = tags$div(sprintf('Whoops. It appears your input PPI file still does not have three columns. Make sure your file only has group/cluster assignment in Column 1, miscellaneous information in Column 2, and gene names as Gene Symbols in Column 3.')),
        closeOnClickOutside = TRUE
      )
      error_ctr <-error_ctr + 1
    }
  } else if (jump_file == 'pathway_db'){
    if (ncol(user_db)!=2){
      sendSweetAlert(
        session=session,
        type = 'warning',
        title='Warning:',
        btn_labels = NA,
        text = tags$div(sprintf('Whoops. It appears your input pathway database does not have two columns. Make sure your file only has the pathway name/description in Column 1, and the genes in the pathway separated by "/" in Column 2.')),
        closeOnClickOutside = TRUE
      )
      error_ctr <-error_ctr + 1
    } else {
      if (nrow(user_db) < 100){
        sendSweetAlert(
          session=session,
          type = 'warning',
          title='Warning:',
          btn_labels = NA,
          text = tags$div(sprintf('Whoops. It appears your input pathway database has less than 100 pathways. Make sure your file only has at least 100 pathways to ensure some decent hits.')),
          closeOnClickOutside = TRUE
        )
        error_ctr <-error_ctr + 1
      } 
    }
  }
  }
  returnValue(error_ctr)
}
