
#################################################
####### Quick Functions For Error Handling ######
#################################################

null_error_handle <- function(object,message){
  if (is.null(object)==TRUE){
    error <- sprintf('<b style="font-size:20px;background-color:tomato;">\t%s</b>', message)
    output$ErrorDisplay <- renderText((HTML(error)))
  }
  validate(
    need(is.null(object)==FALSE, "Please select a data set")
  )
}

typeof_error_handle <- function(object, desired_type,message){
  if (typeof(object)!=desired_type){
    error <- sprintf('<b style="font-size:20px;background-color:tomato;">\t%s</b>', message)
    output$ErrorDisplay <- renderText((HTML(error)))
  }
  validate(
    need(typeof(object)==desired_type, "Please select a data set")
  )
}

equals.equals_error_handle <- function(object1, object2,message){
  if (object1!=object2){
    error <- sprintf('<b style="font-size:20px;background-color:tomato;">\t%s</b>', message)
    output$ErrorDisplay <- renderText((HTML(error)))
  }
  validate(
    need(object1==object2, "Please select a data set")
  )
}