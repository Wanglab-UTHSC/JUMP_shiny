#global.R

#install missing ones
#enable bookmarking

if (!require("pacman")) install.packages("pacman")
pacman::p_load("devtools","shiny","dplyr","shinycssloaders","ggplot2","heatmaply","shinyBS","shinyjs","scales","data.table","RColorBrewer","markdown","utils","shinyWidgets",
               "tidyr","shinydashboard","plotly","gplots","DT","cluster","fresh","shinythemes","limma","clusterProfiler","scatterD3","curl","org.Mm.eg.db", "org.Hs.eg.db",
               "org.Rn.eg.db","glue","spsComps","bookdown")

options(shiny.maxRequestSize = 60*1024^2)

# library(shinydashboard)
# library(plotly)
# library(dplyr)
# library(DT)
# library(heatmaply)
# library(data.table)
# library(RColorBrewer)
# library(markdown)
# library(utils)
# library(shinyWidgets)
# library(tidyr)
# library(cluster)
# library(shinycssloaders)
# library(shinyBS)
# library(MASS)
# library(BiocManager)
# library(shinyjs)
# library(fresh)


# ====================================
# This function convert the input of group information to
# a specific format for TCC calculation.
# ====================================

convert2cl <- function(x, df) {
  which(colnames(df) == x)
}


make_summary_for_tcc_result <- function(df){
  # Make some function for create new column
  sum_gene <- function(x, df){
    sum(df$q.value<=x)
  }
  
  # Set different cut-off
  span <- c(0, seq(0.05, 1, 0.05))
  # Calculate gene count under specific cut-off
  deg_in_cutoff <- sapply(span, sum_gene, df)
  # Calculate total gene count
  total_gene <- nrow(df)
  
  # Create table
  df <- data.frame(
    "Cutoff" = sprintf('%#.2f', span),
    "Under_Count" = deg_in_cutoff,
    "Percentage" = paste(round(deg_in_cutoff / total_gene, 4) * 100, "%")
  )
  df <-
    as_tibble(df) %>% mutate(Between_Count = Under_Count - lag(Under_Count)) %>%
    mutate(Count = paste0(Under_Count, "(+", Between_Count, ")"))
  return(df)
}


actionBttnParams <- list(
  size = "sm",
  color = "default",  # Use a custom color
  style = "bordered",
  class = "action-button",  # Apply the custom class
  block = TRUE
)


radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
    $(document).ready(function() {
      setTimeout(function() {
        $('input', $('#", id, "')).each(function(){
          if(this.getAttribute('value') == '", choice, "') {
            opts = $.extend(", options, ", {html: true});
            $(this.parentElement).tooltip('destroy');
            $(this.parentElement).tooltip(opts);
          }
        })
      }, 500)
    });
  ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}
