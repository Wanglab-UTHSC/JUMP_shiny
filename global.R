#global.R

#Acknowledgement: 
#Some code of this program was adapted from a program related under the MIT license:
#TCC-GUI(https://github.com/swsoyee/TCC-GUI)


if (!require("pacman")) install.packages("pacman")
pacman::p_load("devtools","shiny","dplyr","shinycssloaders","ggplot2","heatmaply","shinyBS","shinyjs","scales","data.table","RColorBrewer","markdown","utils","shinyWidgets",
               "tidyr","shinydashboard","plotly","gplots","DT","cluster","fresh","shinythemes","limma","clusterProfiler","scatterD3","curl","org.Mm.eg.db", "org.Hs.eg.db",
               "org.Rn.eg.db","glue","spsComps","bookdown","R.utils")






#action buttons
actionBttnParams <- list(
  size = "sm",
  color = "default",  # Use a custom color
  style = "bordered",
  class = "action-button",  # Apply the custom class
  block = TRUE
)


#tooltip function on required choice
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
options(shiny.maxRequestSize = 60*1024^2)
