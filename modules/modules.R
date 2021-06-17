reportPCAUI <- function(id, label="PCA plot"){
  ns <- NS(id)
  
  tagList(
    actionButton(ns("genplotpca"), 'Generate plot or update plot', icon("chart-bar"),
                 style="color: #fff; background-color: #006E90; border-color: #006E90"),
    plotOutput(ns("pcaplot")) ## your ggplot should be named output$pcaplot
  )
  
}

reportPCAServer <-function(id){ ## maybe add additional parameters here like KEGG, COG, eggNOG options?
  moduleServer(
    id,
    function(input, output, session) {
      peptide_file <- "peptides.txt"
      metadata_file <- "metadata.csv"
      ## function
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
      # reactive for peptides
      peptides<- reactive({
        peptides  <- read.delim(peptide_file, row.names = 1) %>% as.data.frame() %>% 
          select(starts_with("Intensity."))
      }) 
      
      # reactive for metadata
      metadata<- reactive({
        metadata <- read.delim(metadata_file, header=F, sep=",")
      })
      
      
      # reactive for plot data
      pcaplotdata <- reactive({
        metadata <- metadata()
        cond_options <- table(metadata[,3]) %>% as.data.frame()
        ## filtering for Q50 
        cond_opts <- cond_options$Var1
        cond_count <- cond_options$Freq * 0.5
        
        exp_data <- filter_valids(peptides(),
                                  conditions = cond_opts,
                                  min_count = cond_count,
                                  at_least_one = T)  
        
        norm_pep <- estimateSizeFactorsForMatrix(exp_data) 
        norm_exp <- sweep(as.matrix(exp_data), 2, norm_pep, "/")
        log_data <- data.frame(norm_exp) %>% #dplyr::select(starts_with('Intensity')) %>%
          # mutate_all(., funs(log2(1 + .)))
          mutate(across(everything(), ~{log2(1+.x)}))
        numcond <- length(cond_opts)
        seqcond <- 1:numcond
        colours_to_plot <- lacroix_palette("Pamplemousse", n = numcond, type = "continuous")
        
        colnames(metadata) <- c("Samples", "Experiment", "Condition")
        conditions <- metadata$Condition
        
        pca<- prcomp(t(log_data), center=T, scale=F)
        sampleVals<-data.frame(pca$x)
        exprVals<-data.frame(pca$rotation)
        PoV <- (pca$sdev^2/sum(pca$sdev^2))*100
        
        
        coords<-data.frame(sampleVals, Condition = conditions,
                           samplename = rownames(sampleVals))
        numPCs <- 1:length(PoV)
        
        for (i in 1:length(PoV)) {
          percent <- paste0("(", round(PoV[i],2), "%)")
          percentNoBrack <- paste0(round(PoV[i],2), "%")
          name <- paste0("PC", i, "per")
          name2 <- paste0("PC",i,"per_short")
          assign(name, percent)
          assign(name2, percentNoBrack)
        }
        # now we have all the info we need to plot in "output$plotpca"
        list(coords = coords, colours_to_plot = colours_to_plot, PC1per = PC1per, PC2per = PC2per)
      }) #react close
        
      observeEvent(input$genplotpca,{
        pcaplot <- ggplot(pcaplotdata()[["coords"]], aes(x = PC1, y = PC2)) +
          stat_ellipse(geom = "polygon", alpha=.2, aes(color=Condition, fill=Condition)) +
          geom_point(size=5, aes(colour=Condition, shape=Condition)) + 
          scale_color_manual(values=c(pcaplotdata()[["colours_to_plot"]])) +
          scale_fill_manual(values=c(pcaplotdata()[["colours_to_plot"]])) +
          scale_x_continuous(name= paste0("PC1", " ", pcaplotdata()[["PC1per"]]))+ # labels depend on selected PCs
          scale_y_continuous(name= paste0("PC2", " ", pcaplotdata()[["PC2per"]]))+ 
          theme_bw() +
          theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=12))
        output$pcaplot <- renderPlot({
         pcaplot
        })
      })
  
        
   # your reactive val?
    } # function close
  ) #moduleServer close
}