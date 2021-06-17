library(shiny)
library(tidyverse)
library(LaCroixColoR)
library(DESeq2)
source("modules.R")



ui <- fluidPage(
  reportPCAUI("analysis")
)

server <- function(input, output, session) {
  analysis <- reportPCAServer("analysis")
}

shinyApp(ui, server)