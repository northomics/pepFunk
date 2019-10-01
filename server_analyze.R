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
                label = paste("Input Condition", values$btn + 1), value = ""
      )
    )
  )
  
})

#observe(print(input[['mainDesc1']]))

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
    need(inFile != "", "Please upload a dataset.")
  )
  exp_data <- read.delim(inFile$datapath, row.names = 1) %>% 
    as.data.frame() %>% dplyr::select(starts_with('Intensity.'))
  exp_data[exp_data==1] <-NA

  exp_data
})



output$OriData <-renderRHandsontable({
  if (values$btn > 0) {
    additional_conds <- additional_conds()
    condition_options <- c(input$control, input$othercond, additional_conds(), "NA")
  } else {
    condition_options <- c(input$control, input$othercond, "NA")
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
    
    exp_data <- get_data()
    samplenames <- colnames(exp_data) %>% substr(., 11, nchar(.))
    x <- data.frame(Samples = as.character(samplenames), Condition = as.character(rep("NA", length(samplenames))), 
                    stringsAsFactors = FALSE)
   
     rhandsontable(x) %>%
      hot_col(col = "Condition", type = "dropdown", source = condition_options, strict=T) # must chose a condition
  }
} 
)




#observeEvent(input$runButton, {
#  values$data <-  hot_to_r(input$OriData)
#  
#})

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
style="float:right; color: #fff; background-color: #337ab7; border-color: #2e6da4")
})

get_plotdata <- reactive({
  if (is.null(get_data())) {
    return()
  }
  ## use this for help  
  #  https://deanattali.com/blog/building-shiny-apps-tutorial/ 
  
  exp_data <- get_data()  
  
  ## allow users to ignore samples
  new_conditions <- values$data
  
  if (any("NA"==new_conditions$Condition)) { # only subset the dataframe if use input NA
    ignore_cols <- which(new_conditions$Condition == "NA")
    exp_data <- exp_data[,-ignore_cols]
  }
  
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
    merge(., core_pep_kegg, by='pep', all.x=T) %>% 
    mutate(prop=replace(prop, is.na(prop), 1)) %>%
    mutate_each(funs(.*prop), starts_with('Intensity')) %>% #multiplies the intensities by the proportion 
    mutate(correct_pep = case_when(is.na(newpep_name) ~ pep,
                                   !is.na(newpep_name) ~ newpep_name)) %>%
    column_to_rownames(., var='correct_pep') %>%
    dplyr::select(starts_with('Intensity'))
  
  colnames(core_kegg) <-  removeintensity
  norm_pep <- estimateSizeFactorsForMatrix(core_kegg) 
  exp_data <- sweep(as.matrix(core_kegg), 2, norm_pep, "/")
  peptides <- rownames(exp_data)
  exp_data <- data.frame(exp_data) %>% #dplyr::select(starts_with('Intensity')) %>%
    mutate_all(., funs(log10(1 + .))) # %>% ##should be log10 data...
  rownames(exp_data) <- peptides
  
  
  # applying function over our pathway list
  kegg_genesets <- lapply(pathway_kegg, match_pathway, annot_type='kegg', core_pep_kegg = core_pep_kegg) 
  
  
  ## need to be able to change this
  control_cond <- input$control
  ## make into a function
  
  gsva_kegg <- gsva(as.matrix(exp_data),kegg_genesets, min.sz=10,
                    kcdf='Gaussian') ## rnaseq=F because we have continuous data
  
  
  new_samples <- new_conditions$Samples[-ignore_cols]
  print(new_samples)
  conditions <- new_conditions$Condition[-ignore_cols]
  print(conditions)
  cond <- factor(conditions) %>% relevel(control_cond) # DMSO is the control
  #print(cond)
  design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
  colnames(design)[1] <- c(control_cond) 
  colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                             nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
  fit <- lmFit(gsva_kegg, design)
  fit <- eBayes(fit, trend=T)
  
  colnames(exp_data) <- new_samples
  new_conditions <- data.frame(new_samples, conditions)
  list(exp_data = exp_data, fit = fit, design = design, gsva_kegg = gsva_kegg, conditions = conditions,
       new_conditions = new_conditions)
})

## for the PCA
output$colourpickers <- renderUI({
  # number of conditions will equal to num value$btn + 2
  numcond<- values$btn + 2
  seqcond <- 1:numcond
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
  # number of conditions will equal to num value$btn + 2
  numcond<- values$btn + 2
  seqcond <- 1:numcond
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
    pval <- 1.1 ## change this eventually
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
    melt(., id='Pathway') %>% merge(., new_conditions, by.x='variable', by.y = 'new_samples')
  print(gsvaplot_data %>% head())

  
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
      xlab(label = "Sample") +
      ylab(label="") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) #, plot.margin = margin(6,.8,6,.8, "cm"))
    
  } else {
    values$plotheat <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
      facet_grid(~ conditions, switch='x', scales = "free") +
      #scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
      scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                           space = "Lab", name="GSVA enrichment score") + 
      geom_tile(na.rm = TRUE) +
      xlab(label = "Sample") +
      ylab(label="") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) #, plot.margin = margin(6,.8,6,.8, "cm"))
  }
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
  if (values$btn > 0) {
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
  print(condcolours$Colour %>% as.character())
  dend <- log_exp %>% t() %>% dist(method = dist_method) %>% 
    hclust(method = hclust_method) %>% as.dendrogram() %>%
    set("leaves_pch", 19) %>% 
    set("leaves_col", as.character(condcolours$Colour), order_value = T)
  dend <- as.ggdend(dend, horiz=T)  
  values$dendro <- ggplot(dend)
}) 

## plotting 
output$clustDendro <- renderPlotly({
  validate(
    need(input$genclustdendro, "Please push button to cluster samples and plot or update dendrogram.")
  )
    ggplotly(values$dendro)})
output$heatmapPlot <- renderPlotly({
  
  validate(
    need(input$genplotheat, "Please push button to start analysis and generate or update heatmap.")
  )
  
  ggplotly(values$plotheat, height = 750, width=700})

output$pcaPlot <- renderPlotly({
  validate(
    need(input$genplotpca, "Please push button to start analysis and generate or update PCA biplot.")
  )
  ggplotly(values$plotpca)})

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
    paste0('pcaPlot','.png',sep='')
  }, 
  content = function(file){
    ggsave(file,plot=values$plotpca)
  })

output$downloadKEGG <- downloadHandler(
  filename = function(){
    paste0('peptide_annotation', '.txt')
  },
  content = function(file){
    write.table(core_pep_kegg_only, file, row.names = F, quote = F)}
)