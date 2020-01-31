library(shiny)
fluidPage(
  use_waiter(), # this doesn't show anything, just let's you easily show loading screens
          fluidRow(
                    box(title = "Peptide data input",
                        status = "primary", 
                        width =  4,
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        radioButtons("databaseChoice", "Peptide-to-KEGG database:",
                                     c("Curated human microbiome" = "curated",
                                       "Upload your own database" = "uploaded"),
                                     selected = "curated"),
                       
                        
                        
                        conditionalPanel(
                                         condition = "input.databaseChoice == 'uploaded'", 
                                         tags$hr(),
                                         fileInput("databasefile", "Choose your database file:",
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                       helpText("Note: Your database must be comma separated (.csv) where the first column is the peptide sequence 
                                the second column is the KO term, and the third column is the number of times that peptide was 
                                annotated with the KO term.")
                        ),
                        
                        # horizontal line
                        tags$hr(),
                        
                        radioButtons("sample_data", "Input data type:",
                                     c("Upload your own data" = "upload",
                                       "Use our sample data" = "sample"),
                                     selected = "upload"),
                        
                        # only show if user wants to upload their own data
                        conditionalPanel(
                          condition = "input.sample_data == 'upload'",
                          tags$hr(),
                          radioButtons("file_fmt", "File format:",
                                       c("peptide.txt" = "pep",
                                         "CSV with peptide sequences and intensities" = "csv")),
                          tags$hr(),
                          fileInput("file1", "Choose the peptide intensity file to be analyzed",
                                    accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")
                                    
                          )),
                        
                        #horizontal line
                        tags$hr(),
                        
                       conditionalPanel(
                         condition = "input.sample_data == 'upload'",
                       textInput("control", "Input control condition", placeholder = "Enter control/reference condition name"),
                        textInput("othercond", "Input condition 1", placeholder = "Enter test condition name"),
            
                        tags$div(id = 'placeholder'),
                        actionButton("addcond", "Add additional condition"), 
                        actionButton("rmvcond", "Remove added condition"), # these are ugly
                        tags$hr(),
                        radioButtons("format", "Manual or auto condition formatting?",
                                     c("Manual" = "manual",
                                       "Auto" = "auto"), selected="auto")  ,
                        tags$hr(),
                        radioButtons("normalize", "Transform intensity values using:",
                                     c("Log10" = "log10",
                                       "Log2" = "log2"))
                        #actionButton("runButton","Set condition information")
                ),
                 conditionalPanel(
                   condition = "input.sample_data == 'sample'",
                   textInput("control_sample", "Input control condition", value = "DMSO"),
                   textInput("othercond_sample", "Input condition 1", value = "Low"),
                   textInput("finalcond_sample", "Input condition 2", value = "High")
                 
                 )),
                 
                 
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
                  
