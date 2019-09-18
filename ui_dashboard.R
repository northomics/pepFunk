library(shiny)
fluidPage(
          fluidRow(
                    box(title = "Peptide data input",
                        status = "primary", 
                        width =  4,
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        fileInput("file1", "Choose peptide.txt MaxQuant output",
                                  accept = c(
                                    "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")
                        ),
                        
                        #horizontal line
                        tags$hr(),
                        textInput("control", "Input control condition", placeholder = "Enter control/reference condition name"),
                        textInput("othercond", "Input condition 1", placeholder = "Enter test condition name"),
                        #conditionalPanel(
                        #  condition = "input.moreconditions == 'yes'",
                        #  textInput("othercond2", "Input other condition", "Enter other test condition name")),
                        #radioButtons("moreconditions", "Do you have more conditions?",
                        #             c("Yes" = "yes",
                        #               "No" = "no"), selected="no"),
                        #allow additional conditions, get this working after...
                        tags$div(id = 'placeholder'),
                        actionButton("addcond", "Add additional condition"), 
                        actionButton("rmvcond", "Remove added condition"), # these are ugly
                        tags$hr(),
                        radioButtons("format", "Manual or auto condition formatting?",
                                     c("Manual" = "manual",
                                       "Auto" = "auto"), selected="manual") # ,
                        #tags$hr(),
                        #actionButton("runButton","Set condition information")
                 ),
                 box(title = "Editing sample names and conditions",
                     status = "primary", 
                     width =  8,
                     solidHeader = TRUE,
                     collapsible = FALSE,
                        rHandsontableOutput('OriData')
                 ),
                 
#                 conditionalPanel(
#                   condition = "output.fileUploaded1 == TRUE",
                   box(status = "primary", 
                       width =  8,
                       solidHeader = FALSE,
                       collapsible = FALSE, 
                       uiOutput("gotoanalysisbutton"))
#                       actionButton("gotoanalysis", 
#                                    icon = icon("arrow-right"),
#                                    label = "Continue to peptide centric analysis!",
#                                    style="float:right; color: #fff; background-color: #337ab7; border-color: #2e6da4"
#                       )
                  
            
                 
          )
)
                  
