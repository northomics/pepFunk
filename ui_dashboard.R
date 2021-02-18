library(shiny)
fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
          fluidRow(
         column(5,
                    box(title = "1. Data input",
                        status = "primary", 
                        width = 12,
                        solidHeader = TRUE,
                        collapsible = FALSE,
                        
                       h4( "A. Import peptide file:"), 
                    
                        radioButtons("sample_data", "Input data type:",
                                     c("Upload your own data" = "upload",
                                       "Use our sample data" = "sample"),
                                     selected = "upload"),
                        
                        # only show if user wants to upload their own data
                        conditionalPanel( #begin of conditionalPanel1
                          
                          condition = "input.sample_data == 'upload'",
                          fileInput("file1", "Choose the peptide intensity file to be analyzed",
                                    accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")
                                    
                          ),
                          radioButtons("file_fmt", "File format:",
                                       c("peptide.txt" = "pep",
                                         "CSV with peptide sequences and intensities" = "csv")),
                          conditionalPanel(
                            condition = "input.file_fmt == 'csv'",
                            helpText("File format: A comma separated file (.csv) where the first column is the peptide sequence, and following columns are intensity values.
                                     Headers should be used. It is helpful to have the conditions of your samples in the header.")
                          ),
                          conditionalPanel(
                            condition = "input.file_fmt == 'pep'",
                            helpText("File format: The peptide.txt output file from MaxQuant or MetaLab.")
                          ),
                          ), #end of conditionalPanel1
                       
                       
                       tags$hr(),
                       conditionalPanel(
                         condition = "input.sample_data == 'upload'",
                         h4("B. Add treatment information"),
                         radioButtons("format", "Manual or auto condition formatting?",
                                      c("Manual" = "manual",
                                        "Auto" = "auto"), selected="auto"),
                         helpText("Auto condition formatting will try and match your treatment conditions with your samples. 
                                  If your sample names contain your treatment, this is a great option. 
                                  Try typing your control and treatment names in the boxes above.
                                  If you have more than one treatment, please push the button to add another treatment option."),
                         textInput("control", "Input control condition", placeholder = "Enter control/reference condition name"),
                         textInput("othercond", "Input condition 1", placeholder = "Enter test condition name"),
                         
                         tags$div(id = 'placeholder'),
                         actionButton("addcond", "Add additional condition"), 
                         actionButton("rmvcond", "Remove added condition"), # these are ugly
                         conditionalPanel(
                           condition = "input.sample_data == 'sample'",
                           textInput("control_sample", "Input control condition", value = "DMSO"),
                           textInput("othercond_sample", "Input condition 1", value = "Low"),
                           textInput("finalcond_sample", "Input condition 2", value = "High")
                           
                         ),
                         tags$hr(),
                         
                        
                        
                         
  
                        
    
                        #actionButton("runButton","Set condition information")
                )
                 )),
                 
            column(7,     
                 box(title = "2. Check sample names and sample conditions",
                     status = "primary", 
                     width =  12,
                     solidHeader = TRUE,
                     collapsible = FALSE,
                     # horizontal line
                        rHandsontableOutput('OriData'),
                     
                     tags$hr(),
                     helpText("Note: you can update your sample names here. 
                              Condition names are either auto filled or can be typed in. 
                              Please use the drop down options for conditions.")
                 ),
                 
#                 conditionalPanel(
#                   condition = "output.fileUploaded1 == TRUE",
                   box(title = "3. Analysis options",
                     status = "primary", 
                       width =  12,
                       solidHeader = TRUE,
                       collapsible = FALSE, 
                     h4("A. Data Normalization"),
                     radioButtons("sizefact", "Would you like pepFunk to normalize your data by depth?",
                                  c("Yes" = "yes",
                                    "No" = "no"), 
                                  selected = "yes"),
                     helpText("Note: If you opt for no normalization, we highly recommend you normalize 
                              your data using your own methods before uploading to pepFunk.
                              If you'd like to know more about our normalization technique, please see our manuscript
                              (https://doi.org/10.1093/bioinformatics/btaa289)."),
                     h4("B. Choose log transformation"),
                     radioButtons("normalize", "Transform intensity values using:",
                                  c("Log10" = "log10",
                                    "Log2" = "log2",
                                    "No transformation" = "NoTransf"), 
                                  selected = "log2"),
                     h4("C. Choose peptide annotation type"),
                     ## Annotation type choice
                     radioButtons("databaseType", "Annotation Type:",
                                  c("KEGG" = "kegg",
                                    "COG" = "cog",
                                    "eggNOG" = "eggnog"),
                                  selected = "kegg"),
                     ##KEGG chosen
                     conditionalPanel(
                       condition = "input.databaseChoice == 'kegg'",
                       tags$hr(),
                       radioButtons("databaseChoice", "Peptide-to-KEGG database:",
                                  c("Curated human microbiome" = "curated_kegg",
                                    "Upload your own database" = "uploaded_kegg"),
                                    selected = "curated_kegg"),
                       ##custom panel shows if you want to upload your database
                       conditionalPanel(
                         condition = "input.databaseChoice == 'uploaded_kegg'",
                         tags$hr(),
                         fileInput("databasefile", "Choose your database file:",
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       helpText("Note: Your database must be comma separated (.csv) where the first column is the peptide sequence 
                                the second column is the KO term, and the third column is the number of times that peptide was 
                                annotated with the KO term.")
                       
                     ##COG chosen
                     conditionalPanel(
                       condition = "input.databaseChoice == 'cog'",
                       tags$hr(),
                       radioButtons("databaseChoice", "Peptide-to-COG database:",
                                  c("Curated human microbiome" = "curated_cog",
                                    "Upload your own database" = "uploaded_cog"),
                                    selected = "curated_cog"),
                       ##custom panel shows if you want to upload your database
                       conditionalPanel(
                         condition = "input.databaseChoice == 'uploaded_cog'",
                         tags$hr(),
                         fileInput("databasefile", "Choose your database file:",
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       helpText("Note: Your database must be comma separated (.csv) where the first column is the peptide sequence 
                                the second column is the COG term, and the third column is the number of times that peptide was 
                                annotated with the COG term.")
                       
                     ##eggNOG chosen
                     conditionalPanel(
                       condition = "input.databaseChoice == 'eggnog'",
                       tags$hr(),
                       radioButtons("databaseChoice", "Peptide-to-KEGG database:",
                                  c("Curated human microbiome" = "curated_eggnog",
                                    "Upload your own database" = "uploaded_eggnog"),
                                    selected = "curated_eggnog"),
                       ##custom panel shows if you want to upload your database
                       conditionalPanel(
                         condition = "input.databaseChoice == 'uploaded_eggnog'",
                         tags$hr(),
                         fileInput("databasefile", "Choose your database file:",
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       helpText("Note: Your database must be comma separated (.csv) where the first column is the peptide sequence 
                                the second column is the eggNOG term, and the third column is the number of times that peptide was 
                                annotated with the eggNOG term.")
                       
                     
                     tags$hr(),
                       uiOutput("gotoanalysisbutton"))
#                       actionButton("gotoanalysis", 
#                                    icon = icon("arrow-right"),
#                                    label = "Continue to peptide centric analysis!",
#                                    style="float:right; color: #fff; background-color: #337ab7; border-color: #2e6da4"
#                       )
                  
            
                 
          )
))
                  
