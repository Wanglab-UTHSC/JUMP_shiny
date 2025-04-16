library(shiny)

.libPaths( c( "./lib" , .libPaths() ) )

options(shiny.fullstacktrace=TRUE)

captureStackTraces <- function(expr) {
  withCallingHandlers(
    expr,
    error = function(e) {
      if (is.null(attr(e, "stack.trace", exact = TRUE))) {
        attr(e, "stack.trace") <- sys.calls()
        stop(e)
      }
    }
  )
}

options(shiny.fullstacktrace=TRUE) 
options(shiny.trace = TRUE)
runApp(launch.browser=0, port=3978)
