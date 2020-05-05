#### renderUI 
## allow additional conditions, for adding and removing...
#https://www.reddit.com/r/rstats/comments/7n4qnj/shiny_observeevent_on_inserted_ui/

## loading screen


plotwidth <- reactive({
  plotwidth <- paste0(input$plotwidth)
return(plotwidth)
})
plotheight <- reactive({
  plotheight <- paste0(input$plotheight)
  return(plotheight)
})

## if users want to use their own database
core_pep_kegg <- reactive({
  if (input$databaseChoice == 'curated'){
    ## Core kegg database
    core_pep_kegg <- read.delim("./www/core_pep_kegg_db.csv", 
                                sep=",", header=F, col.names = c("pep", "kegg", "count"))
    core_pep_kegg_only <- core_pep_kegg %>% dplyr::group_by(pep) %>% dplyr::select(pep, kegg)
    core_pep_kegg <- core_pep_kegg %>% dplyr::group_by(pep) %>%
      dplyr::summarize(total = sum(count))  %>%
      merge(., core_pep_kegg, by='pep', all.y=T) %>%
      dplyr::mutate(prop = count/total) %>% dplyr::select(pep, kegg, prop)
    core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep,unique=T) # update pep names so all unique
  } else {
    inFile <- input$databasefile
  validate(
    need(inFile != "", "Please upload a properly formatted peptide to KEGG database.")
  )
  core_pep_kegg <- read.delim(inFile$datapath, 
                              sep=",", header=F, col.names = c("pep", "kegg", "count"))
  core_pep_kegg <- core_pep_kegg %>% dplyr::group_by(pep) %>%
    dplyr::summarize(total = sum(count))  %>%
    merge(., core_pep_kegg, by='pep', all.y=T) %>%
    dplyr::mutate(prop = count/total) %>% dplyr::select(pep, kegg, prop)
  core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep,unique=T) # update pep names so all unique
  }
  return(core_pep_kegg)
})
## need to change in code below
## core_pep_kegg is now core_pep_kegg()

values <- reactiveValues(btn = 0) # want to start the button count at 0 not 1... 


observeEvent(input$addcond, {
  values$btn <- values$btn + 1
  insertUI(
    selector = '#placeholder',
    where = "beforeBegin",
    ui = tags$div( #wrapping in a div with id for ease of removal
      id = paste0('line', values$btn),
      textInput(paste0("otherConditions", values$btn + 1),
                label = paste("Input condition", values$btn + 1), value = ""
      )
    )
  )
  
})


## allow removal of added condition
observeEvent(input$rmvcond, {
  if (values$btn > 0){
  removeUI(
    ## pass in appropriate div id
    selector = paste0('#line', values$btn)
  )
  values$btn <- values$btn - 1
  } else {
    showNotification(
      "You cannot remove any more conditions.",
      duration = 5,
      type = "error") 
  }
})  



additional_conds <- reactive({
  msg <- c(input[["otherConditions"]])
  if (values$btn > 0) {
    for (i in 1:values$btn) {
      msg <- c(as.character(msg), as.character(input[[paste0("otherConditions", i + 1)]]))
    }
    new_conds <- msg
    return(new_conds)
  }
})

get_data <- reactive({
  if (input$sample_data == "sample"){
    inFile <- list(datapath = "peptides.txt")
    } else {
  inFile <- input$file1 }
  validate(
    need(inFile != "", "Please upload a dataset.")
  )
  
#    exp_data <- read.delim("", row.names = 1) %>% 
#      as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
#    exp_data[exp_data==1] <-NA
#  }
  if (input$file_fmt == "pep"){
    exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
      as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
    exp_data[exp_data==1] <-NA
  } else {
    exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
      as.data.frame()
    exp_data[exp_data==1] <-NA
    values$btn <- 0 #resetting btn
  }
  exp_data
})



output$OriData <-renderRHandsontable({
  if (values$btn > 0) {
    additional_conds <- additional_conds()
    condition_options <- c(input$control, input$othercond, additional_conds(), NA)
  } else if (input$sample_data == 'sample') {
    condition_options <-  c(input$control_sample, input$othercond_sample, input$finalcond_sample, NA)
  } else {
    condition_options <- c(input$control, input$othercond, NA)
  }
 
  values$condition_options <- condition_options
  if (input$format == "auto" | input$sample_data == "sample"){
    exp_data <- get_data()
    samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
    x <- data.frame(Samples = samplenames)
    conditions <- purrr::map(condition_options, 
                             ~quo(str_detect(samplenames, fixed(!!.x, ignore_case = T))~!!.x))
    x <- x %>% mutate(Condition = case_when(!!!conditions))
    rhandsontable(x) %>%
      hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
  } else {
    
    exp_data <- get_data()
    samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
    x <- data.frame(Samples = as.character(samplenames), Condition = as.character(rep(NA, length(samplenames))), 
                    stringsAsFactors = FALSE)
   
     rhandsontable(x) %>%
      hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
  }
} 
)


observeEvent(
  input$gotoanalysis, {
    values$data <-  hot_to_r(input$OriData)
    updateTabItems(session, "tabs", "Analysis")
  }
)


output$gotoanalysisbutton <- renderUI({
  if (is.null(get_data())) {
    return()
  }
  actionButton("gotoanalysis", 
icon = icon("arrow-right"),
label = "Continue to peptide centric analysis!",
style="float:right; color: #fff;  background-color: #006E90; border-color: #006E90")
})


output$control_gsva <- renderUI({
  x <- values$data
  conditions <- x$Condition
  if (input$restrict_analysis == "y"){
    selectInput("control_gsva_select", "Control condition",
               choices=conditions)
    }
})

output$treatment_gsva <- renderUI({
  x <- values$data
  conditions <- x$Condition
  if (input$restrict_analysis == "y"){
    selectInput("treatment_gsva_select", "Treatment condition",
                choices=conditions)
  }
})

get_plotdata <- reactive({
  
  withProgress(message = 'Completing GSVA and data transformation', value = 0, { #want a progress bar
  if (is.null(get_data())) {
    return()
  }

  ## use this for help  
  #  https://deanattali.com/blog/building-shiny-apps-tutorial/ 
  core_pep_kegg <- core_pep_kegg()
  exp_data <- get_data()  
  ## allow users to ignore samples
  new_conditions <- values$data
  if (any(is.na(new_conditions$Condition == T))) { # only subset the dataframe if use input NA
    ignore_cols <- is.na(new_conditions$Condition)
    exp_data <- exp_data[,-ignore_cols]
  }

  
  removeintensity <- colnames(exp_data) %>% substr(., 11, nchar(.))
  
 
  core_kegg <- exp_data %>% as.data.frame() %>% 
    rownames_to_column(., var='pep') %>%
    merge(., core_pep_kegg, by='pep', all.x=T) %>% 
    mutate(prop=replace(prop, is.na(prop), 1)) %>%
    mutate_each(funs(.*prop), starts_with('Intensity')) %>% #multiplies the intensities by the proportion 
    mutate(correct_pep = case_when(is.na(newpep_name) ~ pep,
                                   !is.na(newpep_name) ~ newpep_name)) %>%
    column_to_rownames(., var='correct_pep') %>%
    dplyr::select(starts_with('Intensity'))
  
  colnames(core_kegg) <-  removeintensity
  
  if (exists("ignore_cols")){
    new_samples <- new_conditions$Samples[-ignore_cols]
    conditions <- new_conditions$Condition[-ignore_cols]
  } else {
    new_samples <- new_conditions$Samples
    conditions <- new_conditions$Condition
  }
  
  
  cond_opts <- conditions %>% unique()

  cond_count <- table(conditions) %>% as.vector()

  # filtering/removing missing data
  core_kegg <- filter_valids(core_kegg, #we now have 51,395 (from 101,995)
                           #conditions = c('1', '2', '3', '4', '5', '6'),
                           #min_count = c(22, 23, 23, 22, 24, 24), #want peptide to have been identified in at least half of the samples 
                           conditions = cond_opts,
                           min_count = cond_count,
                           at_least_one = TRUE)  #but if it is consistently identified in a condition, keep it

  norm_pep <- estimateSizeFactorsForMatrix(core_kegg) 
  exp_data <- sweep(as.matrix(core_kegg), 2, norm_pep, "/")
  peptides <- rownames(exp_data)
  if (input$normalize == "log10"){
  exp_data <- data.frame(exp_data) %>% #dplyr::select(starts_with('Intensity')) %>%
    mutate_all(., funs(log10(1 + .)))
  } else{
    exp_data <- data.frame(exp_data) %>% #dplyr::select(starts_with('Intensity')) %>%
      mutate_all(., funs(log2(1 + .)))
  }# %>% ##should be log10 data...
  rownames(exp_data) <- peptides
  
  
  # applying function over our pathway list
  kegg_genesets <- lapply(pathway_kegg, match_pathway, annot_type='kegg', core_pep_kegg = core_pep_kegg) 
  
  gsva_kegg <- gsva(as.matrix(exp_data),kegg_genesets, min.sz=10,
                    kcdf='Gaussian') ## rnaseq=F because we have continuous data
  
}) # <- end of withProgress
  list(exp_data = exp_data, gsva_kegg = gsva_kegg, conditions = conditions,
            new_conditions = new_conditions)
})


## render options for GSVA
  output$control_gsva <- renderUI({
    ## render options for GSVA
    condition_options <- get_plotdata()[["conditions"]]
    if (input$restrict_analysis == "y"){
      select_condition <- factor(condition_options)
      #names(select_condition) <- condition_options
      selectInput("control_gsva_select", "Select control condition:",
                  select_condition)
    }else{
      return()
    }
  })
  
  output$treatment_gsva <- renderUI({
    ## render options for GSVA
    condition_options <- get_plotdata()[["conditions"]]
    if (input$restrict_analysis == "y"){
      select_condition <- factor(condition_options)
      #names(select_condition) <- condition_options
      selectInput("treatment_gsva_select", "Select control condition:",
                  select_condition)
    }else{
      return()
    }
  })
  


## for the PCA
output$colourpickers <- renderUI({
  # number of conditions will equal to num value$btn + 2
  if (input$sample_data == "input"){
  numcond<- values$btn + 2
  seqcond <- 1:numcond
  } else {
    cond_options <- values$condition_options
    numcond <- cond_options[!is.na(cond_options)] %>% length()
    values$btn <- 1
    seqcond <- 1:numcond
  }
  # make label for number of conditions
  condition_label <- c("Colour for control/reference", "Colour for condition 1")
  if (values$btn > 0) {
    additional_label <- paste("Colour for condition", 2:values$btn)
    condition_label <- c(condition_label, additional_label)
  }
  ## how to make a palette..
  #colours_to_plot <- viridis_pal(option = "D")(numcond)
  colours_to_plot <- lacroix_palette("Pamplemousse", n = numcond, type = "continuous")
  lapply(seq(numcond), function(i) {
    colourInput(inputId = paste0("colour", seqcond[i]), 
                label = condition_label[i],
                # label = need to figure out
                value = colours_to_plot[i])
     })
})

# for the dendrogram
output$colourpickers2 <- renderUI({
  if (input$sample_data == "input"){
    numcond<- values$btn + 2
    seqcond <- 1:numcond
  } else {
    cond_options <- values$condition_options
    numcond <- cond_options[!is.na(cond_options)] %>% length()
    values$btn <- 1
    seqcond <- 1:numcond
  }
  # make label for number of conditions
  condition_label <- c("Colour for control/reference", "Colour for condition 1")
  if (values$btn > 0) {
    additional_label <- paste("Colour for condition", 2:values$btn)
    condition_label <- c(condition_label, additional_label)
  }
  ## how to make a palette..
  colours_to_plot <- lacroix_palette("Pamplemousse", n = numcond, type = "continuous") 
  #colours_to_plot <- viridis_pal(option = "D")(numcond)
  lapply(seq(numcond), function(i) {
    colourInput(inputId = paste0("colour2_", seqcond[i]), 
                label = condition_label[i],
                # label = need to figure out
                value = colours_to_plot[i])
  })
})

pca_plotdata <- reactive({
  log_exp <- get_plotdata()[['exp_data']]
  new_conditions <- get_plotdata()[['conditions']] # getting the condition data from user's manual input
  pca<- prcomp(t(log_exp), center=T, scale=F)
  sampleVals<-data.frame(pca$x)
  exprVals<-data.frame(pca$rotation)
  PoV <- (pca$sdev^2/sum(pca$sdev^2))*100
  coords<-data.frame(sampleVals, condition = new_conditions,
                     samplename = rownames(sampleVals))
  numPCs <- 1:length(PoV)
  list(coords = coords, PoV = PoV, numPCs=numPCs)
})


## heatmapPlot now is its own separate thing
observeEvent(input$genplotheat,{
  gsva_kegg <- get_plotdata()[['gsva_kegg']]
  new_conditions <- get_plotdata()[['new_conditions']]
  new_samples <- new_conditions$Samples
  
  withProgress(message = 'Making plot', value = 0, { #want a progress bar
  
  if (input$restrict_analysis == "y" ){ # make design matrix for restricted analysis (pairwise comparisons)
    control <- input$control_gsva_select
    treatment <- input$treatment_gsva_select
    new_conditions <- data.frame(new_samples, new_conditions)
    new_conditions <- new_conditions[new_conditions$Condition %in% c(control, treatment),]
    cond <- factor(new_conditions$Condition) %>% relevel(control)
    design <- model.matrix(~ cond)
    colnames(design)[1] <- c(control) 
    colnames(design)[2] <- substr(colnames(design)[2], 5, 
                                  nchar(colnames(design)[2])) #just removing "cond"
    fit <- lmFit(gsva_kegg[,new_conditions$new_samples], design)
    fit<- eBayes(fit, trend=T)

  } else { #plot and analyse ALL the data (no restrictions)
    control_cond <- input$control_sample
    cond <- factor(new_conditions$Condition) %>% relevel(control_cond) # DMSO is the control
    design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
    colnames(design)[1] <- c(control_cond) 
    colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                               nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
    fit <- lmFit(gsva_kegg, design)
    fit <- eBayes(fit, trend=T)

  }
  
  allGeneSets <- topTable(fit, coef=2:ncol(design), number=Inf)
  if (input$plotsig == 'y') {
    pval <- as.numeric(input$pvalthresh)
  } else {
    pval <- Inf # infinity limit for pvalues...no restrictions
  }
  DEgeneSets <- topTable(fit, coef=2:ncol(design), number=Inf,
                         p.value=pval, adjust="BH")
  res <- decideTests(fit, p.value=pval, adjust="BH") ## had to meet adjusted pval
  res <- res %>% as.data.frame()
  
  if (ncol(res) > 2) {
    sig_tests <- res[abs(res[,2:ncol(res)]) %>% rowSums(.) > 0,]
  } else {
    sig_tests <- res[abs(res[,2]) > 0,]
  }


  ## Controlling the type of heatmap. 
  ## Bubble plot options
  ## input$fig_type == "bubble"
  ## input$fig_type == "heatmap"
  
 
    sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]
    ## only looking at significantly altered gene sets.
    if (input$restrict_analysis == "y"){
      control_cond <- input$control_gsva_select 
    } else {
      control_cond <- input$control_sample
    }
    
    if (ncol(sig_tests %>% as.data.frame()) >= 2){
      sigpathways <- sig_tests[abs(sig_tests) %>% rowSums(.) > 0,] %>% as.data.frame() %>%
        rownames_to_column(., var='Pathway') %>% dplyr::select(-control_cond)
    } else {
      sigpathways <- as.data.frame(sig_tests %>% abs())   
      sigpathways <- sigpathways[sigpathways > 0,, drop=F] %>% as.data.frame() %>% rownames_to_column(., var='Pathway')
    }
    
    
    if (input$fig_type == "heatmap"){    
    gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
      melt(., id='Pathway') %>% merge(., new_conditions, by.x='variable', by.y = 'Samples')
    
    
    
    
    ## chosing if we want to plot kegg by p-value or by clustering!
    if (input$kegg_ord == 'clust'){
      kegg_order <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
    } else {
      kegg_order <- allGeneSets[order(-allGeneSets$P.Value),] %>% rownames()
    }
    
    gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = kegg_order)
    # gsvaplot_data <- gsvaplot_data %>% filter(Condition != 'NA')
    
    if (input$sample_ord == 'clust'){
      sample_order <- rownames(sig_gsva %>% t())[hclust(dist(sig_gsva %>% t()))$order]
      gsvaplot_data$variable<- factor(gsvaplot_data$variable, levels = sample_order)
      values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
        #facet_grid(~ Condition, switch='x', scales = "free") +
        #scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
        scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                             space = "Lab", name="GSVA enrichment score") + 
        geom_tile(na.rm = TRUE) +
        xlab(label = "\n\n Sample") +
        ylab(label="") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) #, plot.margin = margin(6,.8,6,.8, "cm"))
      
    } else {
      values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
        facet_grid(~ Condition, switch='x', scales = "free") +
        #scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
        scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                             space = "Lab", name="GSVA enrichment score") + 
        geom_tile(na.rm = TRUE) +
        xlab(label = "\n\n Sample") +
        ylab(label="") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) #, plot.margin = margin(6,.8,6,.8, "cm"))
   
    }
  } else if (input$fig_type == "bubble"){
## bubble heatmap prints, but is not right somewhere...colours change depending on plot type. could be an issue with plotly and fill
  wantedRows <- data.frame(Pathway = rownames(DEgeneSets), pvals = DEgeneSets$adj.P.Val)
  
  #gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
  #    melt(., id='Pathway') %>% merge(., new_conditions, by.x='variable', by.y = 'Samples')  
  sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]
  sig_gsva_plotting <- merge(sig_gsva, wantedRows, by.x=0, by.y="Pathway")
  
  
  gsvaplot_data <- data.frame(sig_gsva_plotting) %>% dplyr::rename(Pathway = Row.names) %>%
    melt(., id=c('Pathway', 'pvals')) %>% merge(., new_conditions, by.x='variable', by.y = 'Samples')
 
  #gsvaplot_data$condition <- substr(gsvaplot_data$variable, 1, nchar(as.character(gsvaplot_data$variable))-1)
  if (input$kegg_ord == 'clust') {
    clusterdata <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
    gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
  } else {
    pvalorder <- allGeneSets[order(-allGeneSets$P.Value),] %>% rownames()
    gsvaplot_data$Pathway <- factor(gsvaplot_data$Pathway, levels = pvalorder)
    }
  if (input$sample_ord == 'clust'){
  values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill=value)) +
      #facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours=c(input$low_col, input$low_col, "white", input$high_col, input$high_col),
                         space = "Lab", name="GSVA enrichment score") + 
      geom_point(na.rm = TRUE, shape=21, colour="darkgrey",  aes(size = abs(value))) +
      xlab(label = "Sample") +
      ylab(label="") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  } else {
    values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill=value)) +
      facet_grid(~ Condition, switch = "x", scales = "free_x", space = "free_x") +
      scale_fill_gradientn(colours=c(input$low_col, input$low_col, "white", input$high_col, input$high_col),
                             space = "Lab", name="GSVA enrichment score") + 
      geom_point(na.rm = TRUE, shape=21,colour="darkgrey",  aes(size = abs(value))) +
      xlab(label = "Sample") +
      ylab(label="") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(size=guide_legend(title="Absolute GSVA enrichment score")# +
    
    ) 
  }
  }
  
      }) ## <- this is end of withProgress
    }) 


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
  
  ## colours are stored under input$colours1 : input$coloursn where n is values$btn + 2
  plotcolours <- paste0("input$colour", 1:(values$btn+2)) 
  ## to change the character vector into object names
  plotcolours <- unlist(lapply(plotcolours,function(s) eval(parse(text=s)))) 
  ## prefer to use shapes that can be filled in...
  if (values$btn < 4) {
    shapes <- c(21,22,23,24,25)
  } else {
    shapes <- c(1:25)
  }
  shapes2use <- shapes[1:(values$btn+2)]
  values$plotpca <- ggplot(coords, aes_string(x = paste0('PC', xaxis), y = paste0('PC', yaxis))) + #accept selectInput to choose axes!
    geom_point(size=3, aes_string(fill="condition", shape="condition")) + 
    stat_ellipse(geom = "polygon", alpha=.2, aes_string(color="condition", fill="condition")) +
    #scale_color_manual(values=c(input$control_col, input$cond1_col, input$cond2_col)) + #pick colours for colour picker
    scale_color_manual(values=plotcolours) +
    scale_fill_manual(values=plotcolours) +
    #scale_shape_manual(values=c(22, 21, 24)) +
    scale_shape_manual(values=shapes2use) +
    scale_x_continuous(name=xaxislabel) + # labels depend on selected PCs
    scale_y_continuous(name=yaxislabel) +
    theme(legend.position = "bottom", legend.title = element_blank()) 
  
}) 

observeEvent(input$genclustdendro, {
  if (input$sample_data == 'sample') {
    condition_options <-  c(input$control_sample, input$othercond_sample, input$finalcond_sample)
    values$btn <- 1
  } else if (values$btn > 0) {
    additional_conds <- additional_conds()
    condition_options <- c(input$control, input$othercond, additional_conds())
  } else {
    condition_options <- c(input$control, input$othercond)
  }
  ## colours are stored under input$colours1 : input$coloursn where n is values$btn + 2
  plotcolours <- paste0("input$colour2_", 1:(values$btn+2)) 
  ## to change the character vector into object names
  plotcolours <- unlist(lapply(plotcolours,function(s) eval(parse(text=s)))) 
  log_exp <- get_plotdata()[['exp_data']]
  dist_method <- input$dist_method
  hclust_method <- input$hclust_method
  dd <- dist(log_exp %>% t(), method = dist_method) # be able to chose distance
  hc <- hclust(dd, method = hclust_method) # be able to chose method
  new_conditions <- values$data # this is a dataframe with Samples, Conditions
  new_samples <- new_conditions$Samples
  condcolours <- data.frame(Condition = condition_options, Colour = plotcolours)
  condcolours <- merge(new_conditions, condcolours, by = "Condition")
  dend <- log_exp %>% t() %>% dist(method = dist_method) %>% 
    hclust(method = hclust_method) %>% as.dendrogram(hang=0.1) %>%
    set("leaves_pch", 19) %>% 
    set("leaves_col", as.character(condcolours$Colour), order_value = T) %>%
    set('branches_lwd', 0.6) %>%
    set('labels_cex', 1)
  dend <- as.ggdend(dend, horiz=T)  
  values$dendro <- ggplot(dend,theme = theme_dendro(), offset_labels = -20) + coord_flip() +
    theme(legend.position="none")
}) 

## plotting 
output$clustDendro <- renderPlotly({
  validate(
    need(input$genclustdendro, "Please push button to cluster samples and plot or update dendrogram.")
  )
   ggplotly(values$dendro)
  })


output$heatmapUI <- renderUI({
  output$heatmapPlot <- renderPlot({
    validate(
      need(input$genplotheat, "Please push button to start analysis and generate/update heatmap.")
    )
  
    values$plotheat})
   plotOutput("heatmapPlot", height = input$plotheight, width=input$plotwidth)
})



output$pcaPlot <- renderPlotly({
  validate(
    need(input$genplotpca, "Please push button to start analysis and generate or update PCA biplot.")
  )
  
  ggplotly(values$plotpca) 
  })

## organizing plot download handlers
output$downloadPlot <- downloadHandler(
  filename = function(){
    paste('heatmapPlot','.pdf',sep='')
  }, 
  content = function(file){
    ggsave(file,plot=values$plotheat, height=input$plotheightsave, width=input$plotwidthsave, units="in")
  }
)


output$dlPCA <- downloadHandler(
  filename = function(){
    paste0('pcaPlot','.png',sep='')
  }, 
  content = function(file){
    ggsave(file,plot=values$plotpca)
  })

#output$downloadKEGG <- downloadHandler(
#  filename = function(){
#    paste0('peptide_annotation', '.txt')
#  },
#  content = function(file){
#    core_pep_kegg <- core_pep_kegg()
#    write.table(core_pep_kegg, file, row.names = F, quote = F)}
#)
output$downloadKEGG <- downloadHandler(
  filename = "peptide_annotation.txt",
  content = function(file){
    #core_pep_kegg <- core_pep_kegg()
    write.table(core_pep_kegg(), file, row.names = F, quote = F)}
)
