library(shiny)

.libPaths( c( "./lib" , .libPaths() ) )

runApp(launch.browser=0, port=3978)
