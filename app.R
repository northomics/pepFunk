library(rhandsontable)
library(shiny)
library(shinyWidgets)
library(colourpicker)
library(reshape2)
library(DT)
library(tidyverse)
library(plyr)
library(DESeq2)
library(GSVA)
library(limma)
library(ggdendro)
library(dendextend)

###########
## TO DO ##
###########
## 1. separate functions from all other code 
## 2. Rename "core_kegg" to something more universal
## 3. make conds/design more universal/customizable
## 4. only an option to find significance...

###########
# Options #
###########
options(shiny.maxRequestSize=125*1024^2)
options(stringsAsFactors = FALSE)

# batch effect correction
#source('remove_batcheffect_modified_functions.R')

# uncomment after debugging
#setwd("~/iclouddrive/Documents/shiny_apps/peptide-centric/")


#    exp_data <- read.delim("peptides.txt", row.names=1) %>%
#exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
#           as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
#    exp_data[exp_data==1] <-NA
########
# DATA # # should add files to this folder
########

## Full KEGG database
kegg_L3 <- read.delim("kegg_L3.txt", sep='\t', header=F, 
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc'),
                      colClasses=c('character','character','character','character')) %>% as.data.frame()
pathways <- kegg_L3$L4_desc %>% unique() 
pathways <- pathways[-c(208:229)] #removing KO that are not in brite/pathway
pathway_kegg <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))

## Core kegg database
core_pep_kegg <- read.delim2("core_pep_kegg.csv", 
                             sep=",", header=F, col.names = c("pep", "kegg", "count", "eval"))
core_pep_kegg <- core_pep_kegg %>% dplyr::group_by(pep) %>% 
  dplyr::summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  dplyr::mutate(prop = count/total) %>% dplyr::select(pep, kegg, prop)
core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep,unique=T) # update pep names so all unique


#############
# FUNCTIONS #
#############



## Remove peptides that are only sometimes identified...
## Data filtering function
## https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-missing-value-imputation/
#filter_valids = function(df, conditions, min_count, at_least_one = TRUE) {
#  # df = data frame containing LOG2 data for filtering and organized by data type
#  # conditions = a character vector dictating the grouping
#  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
#  #     number of valid values for each condition for retention
#  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
#  #     FALSE means min_count must be met across all conditions for retention
#  df[df==0] <- NA
#  all_names <- colnames(df) 
#  cond.names = lapply(conditions, # Group column names by conditions
#                      function(x) grep(x, all_names, value = TRUE, perl = TRUE))
#  cond.filter = sapply(1:length(cond.names), function(i) {
#    df2 = df[cond.names[[i]]]   # Extract columns of interest
#    df2 = as.matrix(df2)   # Cast as matrix for the following command
#    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
#    sums >= min_count[i]   # Calculates whether min_count requirement is met
#  })
#  if (at_least_one) {
#    df$KEEP = apply(cond.filter, 1, any)
#  } else {
#    df$KEEP = apply(cond.filter, 1, all)
#  }
#  #return(df) # No rows are omitted, filter rules are listed in the KEEP column
#  df[is.na(df)] <- 0 
#  return(df %>%  rownames_to_column(., var='peptides') %>% filter(KEEP) %>% dplyr::select(-KEEP) %>%
#           column_to_rownames(., var='peptides'))  # only keeping rows that meet the criteria!
#}

expression_data <- function(annotationtype, expr_df){
  if (annotationtype=='COG'){
    norm_pep <- estimateSizeFactorsForMatrix(expr_df)
    exp_data <- sweep(expr_df, 2, norm_pep, "/") ##peptides are normalized
    exp_data <- core_drug_cog %>% #dplyr::select(starts_with('Intensity')) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), starts_with('Intensity')) %>% ##should be log10 data...
      column_to_rownames(., var='peptide') 
    return(expr_df)
  } else if (annotationtype=='kegg'){
    norm_pep <- estimateSizeFactorsForMatrix(expr_df)
    exp_data <- sweep(expr_df, 2, norm_pep, "/") ##peptides are normalized
    exp_data <- core_kegg %>% #dplyr::select(starts_with('Intensity')) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), starts_with('Intensity')) %>% ##should be log10 data...
      column_to_rownames(., var='peptide')
    return(expr_df)    
  } else {
    stop("Not an accepted functional annotation type.")
  }}

# finish writing this function
## write a function that will do this for you!  
match_pathway <- function(df, annot_type){
  if (annot_type == 'kegg'){
    subset <- core_pep_kegg[core_pep_kegg$kegg %in% df$L3,]
    subset <- subset$newpep_name
    return(subset)
  } else if (annot_type=='cog'){
    subset <- core_pep_cog[core_pep_cog$cog %in% df$cog,]
    subset <- subset$newpep_name
    return(subset)
  } else {
    stop("Not an accepted functional annotation type.")}
}

# Function that produces default gg-colours is taken from this discussion:
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## tips
## source scripts, load libraries, and read datasets at the beginning of the script
## outside of the server function...
## shiny will only read this once (which is all you need)

## define user specific objects inside server function but outside of any render calls

## place all code shiny must rerun inside of the render function.
## shiny will rerun all the code in a render chunk each time a user changes a widget...
ui <- navbarPage("Peptide-centric metaproteomic workflow",
                
                  tabPanel("Data setup",
                           fluidRow(
                             column(4, h3("Peptide data input"),
                                    
                                    fileInput("file1", "Choose peptide.txt MaxQuant output",
                                              accept = c(
                                                "text/csv",
                                                "text/comma-separated-values,text/plain",
                                                ".csv"),
                                    ),
                                    
                                    #horizontal line
                                    tags$hr(),
                                    textInput("control", "Input control condition", "Enter control/reference condition name"),
                                    textInput("othercond", "Input other condition", "Enter test condition name"),
                                    conditionalPanel(
                                      condition = "input.moreconditions == 'yes'",
                                      textInput("othercond2", "Input other condition", "Enter other test condition name")),
                                    radioButtons("moreconditions", "Do you have more conditions?",
                                                 c("Yes" = "yes",
                                                   "No" = "no"), selected="no"),
                                    #allow additional conditions, get this working after...
                                    actionButton("addcond", "Add additional condition"), actionButton("rmvcond", "Remove added condition"), # these are ugly
                                    tags$div(id = 'placeholder'),
                                    verbatimTextOutput("value", placeholder = T),
                                    
                                    radioButtons("format", "Manual or auto condition formatting?",
                                                 c("Manual" = "manual",
                                                   "Auto" = "auto"), selected="manual"),
                                    tags$hr(),
                                    actionButton("runButton","Set condition information")
                                    ),
                             column(8, h3("Editing sample names and conditions"),
                                             rHandsontableOutput('OriData'))
                           )),
                 
                 tabPanel("PCA",
                          fluidRow(
                            column(12, plotOutput("pcaPlot")))
                            ,
                          #uiOutput("y_axisPC"), # get back to this but see if this is the issue first
                          fluidRow(
                            column(6,
                          selectInput("y_axisPC", label = "PC on y-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                                      choices = c('2' = '2'),
                                      selected = '2'),
                          selectInput("x_axisPC", label = "PC on x-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                                      choices = c('1' = '1'),
                                      selected = "1"),
                          actionButton('genplotpca', 'Generate PCA biplot and update colours'),
                          downloadButton('dlPCA', 'Download PCA biplot')),
                          column(6,
                                 colourInput("control_col", "Colour for control/reference condition", "#67A7C1"),
                                 colourInput("cond1_col", "Colour for condition 1", "#FF6F59"),
                                 colourInput("cond2_col", "Colour for condition 2", "#292F36")
                                 )
                          
                          )),
                 
                 tabPanel("Sample clustering",
                          
                          fluidRow(
                            column(12,
                          plotOutput("clustDendro")
                            )
                          ),
                          
                          fluidRow(
                            column(6,
                                   colourInput("control_coldend", "Colour for control/reference condition", "#67A7C1"),
                                   colourInput("cond1_coldend", "Colour for condition 1", "#FF6F59"),
                                   colourInput("cond2_coldend", "Colour for condition 2", "#292F36")
                            ),
                            column(6,
                                   selectInput("dist_method", label = "Distance method:",
                                               choices = c('Euclidean' = 'eucl',
                                                        'Canberra' = 'canb',
                                                         'Binary' = 'bina',
                                                        'Minkowski' = 'mink'),
                                                       selected = 'eucl'),
                                   selectInput("hclust_method", label = "Hierarchical clustering method:",
                                               choices = c('Ward' = 'ward.D',
                                                           'Ward 2' = 'ward.D2',
                                                           'Single' = 'sing',
                                                           'Complete' = 'compl',
                                                           'Average' = 'aver',
                                                           'McQuitty' = 'mcq',
                                                           'Median' = 'medi',
                                                           'Centroid' = 'centr'),
                                               selected = 'ward.D2'),
                                   actionButton('genclustdendro', 'Generate cluster dendrogram'),
                                   downloadButton('dlDendro', 'Download cluster dendrogram')
                          )
                          )), #maybe have a second drop down menu for this
                 
                 tabPanel("Functional enrichment heatmap",
                          sidebarLayout(
                            sidebarPanel(width = 3,
                                         colourInput("high_col", "Colour for high GSVA score", "FF6F59"),
                                         colourInput("low_col", "Colour for low GSVA score", "#67A7C1"),
                                         radioButtons("sample_ord", "Order Samples:",
                                                      c('By clustering' = 'clust',
                                                        'By conditions/facets' = 'facets')),
                                         radioButtons("kegg_ord", "Order KEGG pathways:",
                                                      c('By clustering' = 'clust',
                                                        'By p-value' = 'pval')),
                                         radioButtons('plotsig', "Only plot significantly enriched KEGG?",
                                                      c('Yes' = 'y',
                                                        'No' = 'n'),
                                                      selected = 'n'),
                                         conditionalPanel(
                                           condition = "input.plotsig == 'y'",
                                           textInput("pvalthresh", "P-value threshold for plotting", "0.05")),
                                         actionButton('genplotheat', 'Generate/update heatmap'),
                                         downloadButton('downloadPlot','Download heatmap')
                            ),
                            mainPanel(
                              plotOutput("heatmapPlot")
                            ))

                    
))
server <- function(input,output,session)({
#### renderUI 
  ## allow additional conditions, for adding and removing...
  #https://www.reddit.com/r/rstats/comments/7n4qnj/shiny_observeevent_on_inserted_ui/
  
  values <- reactiveValues(btn = 0) # want to start the button count at 0... 
  ## be able to add conditions....testing this!!
  ## this is not working...
#  observeEvent(input$addcond, {
#    condition_options <- c(input$control, input$othercond) 
#    newcond <- as.character(input$add)
#    print(newcond)
#    condition_options <<- c(condition_options, newcond)
#  })
#  observeEvent(input$addcond, {
#    insertUI(
#      selector = "#addcond",
#      where = "beforeBegin",
#      ui = textInput(paste0("txt", input$add),
#                     "Additional condition")
#    )
#  })
  
  observeEvent(input$addcond, {
    values$btn <- values$btn + 1
    insertUI(
      selector = '#placeholder',
      where = "beforeBegin",
      ui = tags$div( #wrapping in a div with id for ease of removal
        id = paste0('line', values$btn),
        textInput(paste0("otherConditions", values$btn + 1),
                  label = paste("Condition", values$btn + 1), value = ""
                  )
       )
    )

  })
  
 #observe(print(input[['mainDesc1']]))
 
  ## allow removal of added condition
  observeEvent(input$rmvcond, {
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#line', values$btn)
    )
    values$btn <- values$btn - 1 
  })  
  
  
  #output$value <- renderText({ 
  additional_conds <- reactive({
    msg <- c(input[["otherConditions"]])
    if (values$btn > 0) {
      for (i in 1:values$btn) {
        msg <- c(as.character(msg), as.character(input[[paste0("otherConditions", i + 1)]]))
      }
      #new_conds <- paste(msg, collapse = '" , "')  
      new_conds <- msg
      return(new_conds)
    }
  })
  
  get_data <- reactive({
    inFile <- input$file1
    validate(
      need(inFile != "", "Please upload a dataset")
    )
    exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
      as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
    exp_data[exp_data==1] <-NA
    exp_data
  })
  
  
  
output$OriData <-renderRHandsontable({
 if (values$btn > 0) {
   additional_conds <- additional_conds()
   condition_options <- c(input$control, input$othercond, additional_conds())
 } else {
   condition_options <- c(input$control, input$othercond)
 }
  
  if (input$format == "auto"){
    exp_data <- get_data()
    samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
    x <- data.frame(Samples = samplenames)
    conditions <- purrr::map(condition_options, 
                             ~quo(str_detect(samplenames, fixed(!!.x, ignore_case = T))~!!.x))
    x <- x %>% mutate(Condition = case_when(!!!conditions))
     rhandsontable(x) %>%
       hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
   } else {
     
     # rhandsontable(as.data.frame(iris))
     exp_data <- get_data()
     #condition_options <- c(input$control, input$othercond, "NA")
     #condition_options <- c("High", "DMSO", "NA")
     samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
     x <- data.frame(Samples = as.character(samplenames), Condition = as.character(rep("NA", length(samplenames))), 
                     stringsAsFactors = FALSE)
     rhandsontable(x) %>%
       hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
   }
      } 
)
  
  
## This works!! I just want to try enabling for "unlimited" conditions...
#  output$OriData <- renderRHandsontable({
#    additional_conds <- additional_conds()
#    if (input$moreconditions=='yes') {
#      condition_options <- c(input$control, input$othercond, input$othercond2, "NA")
#    } else {
#      condition_options <- c(input$control, input$othercond, "NA") 
#    }
#    if (input$format == "auto"){
#      exp_data <- get_data()
#      samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
#      x <- data.frame(Samples = samplenames) %>%
#        mutate(Condition = case_when(
#          str_detect(Samples, fixed(input$control, ignore_case = T)) ~ input$control,
#          str_detect(Samples, fixed(input$othercond, ignore_case = T)) ~ input$othercond,
#          str_detect(Samples, fixed(input$othercond2, ignore_case = T)) ~ input$othercond2,
#          !(str_detect(Samples, fixed(input$control, ignore_case = T)) | str_detect(Samples, fixed(input$othercond, ignore_case = T)))~ "NA"))
#      rhandsontable(x) %>%
#        hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
#    } else {
#      
#      # rhandsontable(as.data.frame(iris))
#      exp_data <- get_data()
#      #condition_options <- c(input$control, input$othercond, "NA")
#      #condition_options <- c("High", "DMSO", "NA")
#      samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
#      x <- data.frame(Samples = as.character(samplenames), Condition = as.character(rep("NA", length(samplenames))), 
#                      stringsAsFactors = FALSE)
#      rhandsontable(x) %>%
#        hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
#    }
#  })
  
  
  observeEvent(input$runButton, {
    values$data <-  hot_to_r(input$OriData)
  })
  
  
  ## Maybe we should make a single reactive for "normal" data manipulation    
  #  pca_data <- reactive({
  #    if (is.null(get_data())) {
  #      return()
  #    }
  #    exp_data <- get_data()
  #    })
  
  
  get_plotdata <- reactive({
    if (is.null(get_data())) {
      return()
    }
    ## use this for help  
    #  https://deanattali.com/blog/building-shiny-apps-tutorial/ 
    
    exp_data <- get_data()  
    removeintensity <- colnames(exp_data) %>% substr(., 11, nchar(.))
    ### will need to uncomment 
    #    # exp_data = filter_valids(exp_data,
    #    #   conditions = c('DMSO', 'High', 'Low'),
    #    #   min_count = c(2, 3, 2), #want peptide to have been identified in at least half of the samples 
    #    #   at_least_one = TRUE)
    #    
    #    
    #    ## intensities normalized by the proportion of functional annotation
    #    ## if there are more than one functional annotation, the peptide will have a suffix added to the end (i.e. .1, .2, .3...etc)
    #    ## $newpep_name
    core_kegg <- exp_data %>% as.data.frame() %>% 
      rownames_to_column(., var='pep') %>%
      merge(., core_pep_kegg, by='pep') %>% 
      mutate(prop=replace(prop, is.na(prop), 1)) %>%
      mutate_each(funs(.*prop), starts_with('Intensity')) %>% #multiplies the intensities by the proportion 
      column_to_rownames(., var='newpep_name') %>%
      dplyr::select(starts_with('Intensity'))
    
    colnames(core_kegg) <-  removeintensity
    norm_pep <- estimateSizeFactorsForMatrix(core_kegg) 
    exp_data <- sweep(as.matrix(core_kegg), 2, norm_pep, "/")
    peptides <- rownames(exp_data)
    exp_data <- core_kegg %>% #dplyr::select(starts_with('Intensity')) %>%
      mutate_all(., funs(log(1 + .))) # %>% ##should be log10 data...
    rownames(exp_data) <- peptides
    
    
    # applying function over our pathway list
    kegg_genesets <- lapply(pathway_kegg, match_pathway, annot_type='kegg') 
    ## shiny app has a UI, server function, then call to the shiny app...
    
    ## need to be able to change this
    control_cond <- input$control
    ## make into a function
    gsva_kegg <- gsva(as.matrix(exp_data),kegg_genesets, min.sz=10,
                      kcdf='Gaussian') ## rnaseq=F because we have continuous data
    new_conditions <- values$data
    new_samples <- new_conditions$Samples
    cond <- factor(new_conditions$Condition) %>% relevel(control_cond) # DMSO is the control
    #print(cond)
    design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
    colnames(design)[1] <- c(control_cond) 
    colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                               nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
    fit <- lmFit(gsva_kegg, design)
    fit <- eBayes(fit, trend=T)
 
    colnames(exp_data) <- new_samples
    list(exp_data = exp_data, fit = fit, design = design, gsva_kegg = gsva_kegg, new_conditions = new_conditions)
  })
  
  pca_plotdata <- reactive({
    log_exp <- get_plotdata()[['exp_data']]
    new_conditions <- values$data # getting the condition data from user's manual input
    pca<- prcomp(t(log_exp), center=T, scale=F)
    sampleVals<-data.frame(pca$x)
    exprVals<-data.frame(pca$rotation)
    PoV <- (pca$sdev^2/sum(pca$sdev^2))*100
    
    # Make is so that it loops through all possibilities
    #for (i in length(PoV)){
    #  
    #}
    
    # need to make this more general
    coords<-data.frame(sampleVals, condition = new_conditions,
                       samplename = rownames(sampleVals))
    #values$numPCs <- paste0("PC", 1:length(PoV))
    numPCs <- 1:length(PoV)
    #print(numPCs)
    ## dropdown for selecting which PC we want to plot
    
    list(coords = coords, PoV = PoV, numPCs=numPCs)
  })
  

    
  ## heatmapPlot now is its own separate thing, want to be able to push a button for this...
  observeEvent(input$genplotheat,{
    fit <- get_plotdata()[['fit']]
    design <- get_plotdata()[['design']]
    gsva_kegg <- get_plotdata()[['gsva_kegg']]
    new_conditions <- get_plotdata()[['new_conditions']]
    allGeneSets <- topTable(fit, coef=2:ncol(design), number=Inf)
    ## only plot "significant" pvalues
    if (input$plotsig == 'y') {
      pval <- as.numeric(input$pvalthresh)
    } else {
      pval <- 0.05 ## change this eventually
    }
    DEgeneSets <- topTable(fit, coef=2:ncol(design), number=Inf,
                           p.value=pval, adjust="BH")
    res <- decideTests(fit, p.value=pval)
    
    sigdrugs <- res[,abs(res) %>% colSums(.) > 0]
    
    control_cond <- input$control
    
    ## only looking at significantly altered gene sets.
    if (ncol(sigdrugs %>% as.data.frame()) >= 2){
      sigpathways <- sigdrugs[abs(sigdrugs) %>% rowSums(.) > 0,] %>% as.data.frame() %>%
        rownames_to_column(., var='Pathway') %>% dplyr::select(-control_cond)
    } else {
      sigpathways <- as.data.frame(sigdrugs %>% abs())   
      sigpathways <- sigpathways[sigpathways > 0,, drop=F] %>% as.data.frame() %>% rownames_to_column(., var='Pathway')
    }
    
    
    sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% sigpathways$Pathway,]
    gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
      melt(., id='Pathway') %>%
      merge(., new_conditions, by.x='variable', by.y = 'Samples')
    
    #clusterdata <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
    
    ## chosing if we want to plot kegg by p-value or by clustering!
    if (input$kegg_ord == 'clust'){
      kegg_order <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
    } else {
      kegg_order <- allGeneSets[order(-allGeneSets$P.Value),] %>% rownames()
    }
    
    gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = kegg_order)
    gsvaplot_data <- gsvaplot_data %>% filter(Condition != 'NA')
    
    if (input$sample_ord == 'clust'){
      sample_order <- rownames(sig_gsva %>% t())[hclust(dist(sig_gsva %>% t()))$order]
      gsvaplot_data$variable<- factor(gsvaplot_data$variable, levels = sample_order)
      values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
        #facet_grid(~ Condition, switch='x', scales = "free") +
        #scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
        scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                             space = "Lab", name="GSVA enrichment score") + 
        geom_tile(na.rm = TRUE) +
        xlab(label = "Sample") +
        ylab(label="") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
    } else {
      values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
        facet_grid(~ Condition, switch='x', scales = "free") +
        #scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
        scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                             space = "Lab", name="GSVA enrichment score") + 
        geom_tile(na.rm = TRUE) +
        xlab(label = "Sample") +
        ylab(label="") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }) 
  
 # yaxis_pcObs <- reactive({
 #   if (is.null(input$y_axisPC)){
 #     yaxis <- "PC2" # how do you make this work?
 #   } else {
 #   yaxis <- input$y_axisPC 
 #   }
 #   print(yaxis)
 # })
  
  
  
  observeEvent(input$genplotpca, {
    coords <- pca_plotdata()[['coords']]
    PoV <- pca_plotdata()[['PoV']]
    numPCs <- pca_plotdata()[['numPCs']]
    for (i in 1:length(PoV)) {
      percent <- paste0("(", round(PoV[i],2), "%)")
      name <- paste0("PC", i, "per")
      assign(name, percent)
      }
    yaxis <- input$y_axisPC
    xaxis <- input$x_axisPC
    updateSelectInput(session, "y_axisPC", label = "PC on y-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                      choices = as.list(numPCs),
                      selected = yaxis)
    updateSelectInput(session, "x_axisPC", label = "PC on x-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                      choices = as.list(numPCs),
                      selected = xaxis)
    yaxis <- input$y_axisPC
    xaxis <- input$x_axisPC
    yperc <- paste0("(", round(PoV[yaxis %>% as.numeric()] ,2), "%)")
    xperc <- paste0("(", round(PoV[xaxis %>% as.numeric()] ,2), "%)")
    
    yaxislabel <- paste0("PC", yaxis, " ", yperc)
    xaxislabel <- paste0("PC", xaxis, " ", xperc)
    values$plotpca <- ggplot(coords, aes_string(x = paste0('PC', xaxis), y = paste0('PC', yaxis))) + #accept selectInput to choose axes!
      geom_point(size=3, aes_string(fill="condition.Condition", shape="condition.Condition")) + 
      stat_ellipse(geom = "polygon", alpha=.2, aes_string(color="condition.Condition", fill="condition.Condition")) +
      scale_color_manual(values=c(input$control_col, input$cond1_col, input$cond2_col)) + #pick colours for colour picker
      scale_fill_manual(values=c(input$control_col, input$cond1_col, input$cond2_col)) +
      scale_shape_manual(values=c(22, 21, 24)) +
      scale_x_continuous(name=xaxislabel) + # labels depend on selected PCs
      scale_y_continuous(name=yaxislabel) +
      theme(legend.position = "bottom", legend.title = element_blank()) 
    
  }) 
  
  observeEvent(input$genclustdendro, {
  #observe({
    log_exp <- get_plotdata()[['exp_data']]
    dist_method <- input$dist_method
    hclust_method <- input$hclust_method
    dd <- dist(log_exp %>% t(), method = dist_method) # be able to chose distance
    hc <- hclust(dd, method = hclust_method) # be able to chose method
    condcolours <- values$data %>% mutate(Colour = case_when(
      Condition == input$control ~ input$control_coldend,
      Condition == input$othercond ~ input$cond1_coldend,
      Condition == input$othercond2 ~ input$cond2_coldend)
      )
    dend <- log_exp %>% t() %>% dist(method = dist_method) %>% 
      hclust(method = hclust_method) %>% as.dendrogram() %>%
      set("leaves_pch", 19) %>%
      set("leaves_col", condcolours$Colour, order_value = TRUE)
    dend <- as.ggdend(dend)  
    values$dendro <- ggplot(dend)
  }) 
  
  ## plotting 
  output$clustDendro <- renderPlot({values$dendro}, height='auto')
  output$heatmapPlot <- renderPlot({values$plotheat}, height='auto')
  output$pcaPlot <- renderPlot({values$plotpca}, height='auto')
  
  ## organizing plot download handlers
  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste('heatmapPlot','.png',sep='')
    }, 
    content = function(file){
      ggsave(file,plot=values$plotheat)
    }
  )
  
  
  output$dlPCA <- downloadHandler(
    filename = function(){
      paste('pcaPlot','.png',sep='')
    }, 
    content = function(file){
      ggsave(file,plot=values$plotpca)
    })
  
})



shinyApp(ui, server)