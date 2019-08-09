library(shiny)
library(tidyverse)
library(plyr)
library(DESeq2)
library(GSVA)
library(limma)
library(reshape2)
###########
## TO DO ##
###########
## 1. separate functions from all other code 
## 2. Rename "core_drug_kegg" to something more universal
## 3. make conds/design more universal/customizable
## 4. only an option to find significance...

###########
# Options #
###########
options(shiny.maxRequestSize=100*1024^2)

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
filter_valids = function(df, conditions, min_count, at_least_one = TRUE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  df[df==0] <- NA
  all_names <- colnames(df) 
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, all_names, value = TRUE, perl = TRUE))
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  #return(df) # No rows are omitted, filter rules are listed in the KEEP column
  df[is.na(df)] <- 0 
  return(df %>%  rownames_to_column(., var='peptides') %>% filter(KEEP) %>% dplyr::select(-KEEP) %>%
           column_to_rownames(., var='peptides'))  # only keeping rows that meet the criteria!
}

expression_data <- function(annotationtype){
  if (annotationtype=='COG'){
    norm_pep <- estimateSizeFactorsForMatrix(exp_data)
    exp_data <- sweep(exp_data, 2, norm_pep, "/") ##peptides are normalized
    exp_data <- core_drug_cog %>% #dplyr::select(starts_with('Intensity')) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), starts_with('Intensity')) %>% ##should be log10 data...
      column_to_rownames(., var='peptide') 
    return(exp_data)
  } else if (annotationtype=='kegg'){
    norm_pep <- estimateSizeFactorsForMatrix(exp_data)
    exp_data <- sweep(exp_data, 2, norm_pep, "/") ##peptides are normalized
    exp_data <- core_drug_kegg %>% #dplyr::select(starts_with('Intensity')) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), starts_with('Intensity')) %>% ##should be log10 data...
      column_to_rownames(., var='peptide')
    return(exp_data)    
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



## tips
## source scripts, load libraries, and read datasets at the beginning of the script
## outside of the server function...
## shiny will only read this once (which is all you need)

## define user specific objects inside server function but outside of any render calls

## place all code shiny must rerun inside of the render function.
## shiny will rerun all the code in a render chunk each time a user changes a widget...

######
# UI #
#######

ui <- fluidPage(
  titlePanel("Peptide-centric functional analysis"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Look for KEGG functional enrichment in peptide intensities
               using Gene Set Varation Analysis (GSVA)."),

      
      fileInput("file1", "Choose CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      
      ## horizontal line
      tags$hr(),
      
     
    
    # Input: Select number of rows to display ----
    radioButtons("disp", "Display",
                 choices = c(Head = "head",
                             All = "all"),
                 selected = "head"),
    
    tags$hr(),
    
    textInput("control_cond", h3("Control Condition"), 
              value = "Name of control condition for plot"),
    
    textInput("other_cond", h3("Other Condition"), 
              value = "Name of other condition for plot")
    
  ), 
   
    mainPanel(
      tableOutput("contents"),
      plotOutput("heatmapPlot")
#      textOutput("file_name")
#      textOutput("min_max")
    )
)
)


# server function is what will actually happen with the data

##########
# SERVER #
##########

server <- function(input, output, session) {
  
  
  get_data<-reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL)
    } # if no upload
    exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
           as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
    exp_data[exp_data==1] <-NA
    return(exp_data)
})

    
  output$contents <- renderTable({
    exp_data <- get_data()
    if(input$disp == "head") {
      return(head(exp_data))
    }
    else {
      return(exp_data)
  }
})
  
#  output$heatmapPlot <- renderPlot({
#    # need to make this so the user can select the column
#    ## need to make this more custom...
#    ## Apply filtering
cond <- colnames(exp_data) %>% substr(., 1, nchar(.)-1)
  #    exp_data = filter_valids(exp_data(),
#                             conditions = c('DMSO', 'High', 'Low'),
#                             min_count = c(2, 3, 2), #want peptide to have been identified in at least half of the samples 
#                             at_least_one = TRUE)
#    
#    
#    ## intensities normalized by the proportion of functional annotation
#    ## if there are more than one functional annotation, the peptide will have a suffix added to the end (i.e. .1, .2, .3...etc)
#    ## $newpep_name
#    core_drug_kegg <- exp_data %>% 
#      rownames_to_column(., var='pep') %>%
#      merge(., core_pep_kegg, by='pep') %>% 
#      mutate(prop=replace(prop, is.na(prop), 1)) %>%
#      mutate_each(funs(.*prop), starts_with('Intensity')) %>% #multiplies the intensities by the proportion 
#      column_to_rownames(., var='newpep_name') %>%
#      dplyr::select(starts_with('Intensity'))
#    
#    
#    ## logged and normalized data
#    exp_data_kegg <- expression_data('kegg')
#  
#    
#    # applying function over our pathway list
#    kegg_genesets <- lapply(pathway_kegg, match_pathway, annot_type='kegg') 
#    ## shiny app has a UI, server function, then call to the shiny app...
#    
#    ## need to be able to change this
#    control_cond <- "DMSO"
#    ## make into a function
#    gsva_kegg <- gsva(as.matrix(exp_data_kegg),kegg_genesets, min.sz=10,
#                      kcdf='Gaussian') ## rnaseq=F because we have continuous data
#    cond<- factor(cond) %>% relevel(control_cond) # DMSO is the control
#    design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
#    colnames(design)[1] <- c(control_cond) 
#    colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
#                                               nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
#    fit <- lmFit(gsva_kegg, design)
#    fit <- eBayes(fit, trend=T)
#    allGeneSets <- topTable(fit, coef=2:ncol(design), number=Inf)
#    DEgeneSets <- topTable(fit, coef=2:ncol(design), number=Inf,
#                           p.value=0.05, adjust="BH")
#    res <- decideTests(fit, p.value=0.05)
#    sigdrugs <- res[,abs(res) %>% colSums(.) > 0]
#    ## hierarcical clustering of the drugs by kegg
#    clusterdata <- colnames(gsva_kegg)[hclust(dist(gsva_kegg%>%t()))$order]
#    
#    ## only looking at significantly altered gene sets.
#    #sigpathways <- sigdrugs[abs(sigdrugs) %>% rowSums(.) > 0,] %>% as.data.frame() %>%
#    #  rownames_to_column(., var='Pathway') %>% dplyr::select(-control_cond)
#    #sigpathways <- melt(sigpathways, id.vars='Pathway') %>% filter(value != 0) %>% dplyr::select(-value) %>%
#    #  mutate(keep = rep("KEEP", nrow(.)))
#    
#    sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(DEgeneSets),]
#    gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
#      melt(., id='Pathway') 
#    gsvaplot_data$condition <- substr(gsvaplot_data$variable, 1, nchar(as.character(gsvaplot_data$variable))-1)
#    
#    clusterdata <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
#    gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
#    
#    (highplot <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
#       facet_grid(~ condition, switch = "x", scales = "free_x", space = "free_x") +
#       scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
#                            space = "Lab", name="GSVA enrichment score") + 
#       geom_tile(na.rm = TRUE) +
#       xlab(label = "Sample") +
#       ylab(label="") +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)))
#  })
  
  #output$file_name <- renderText({ 
  #  paste("You have selected", input$file)
  #})
  
  #output$plot
  
}
# this is the call to the shinyApp
shinyApp(ui, server)
#you would run this app outside the directly as runApp("directoryname")