observeEvent({ 
  input$spec_button
  mainplot.click$click
},
{
  if (input$restrict_analysis == "y"){ # make design matrix for restricted analysis (pairwise comparisons)
    control <- input$control_gsva_select
    print(control)
    treatment <- input$treatment_gsva_select
    print(treatment)
    new_conditions <- data.frame(new_samples, conditions)
    print(new_conditions)
    new_conditions <- new_conditions[new_conditions$conditions %in% c(control, treatment),]
    print(new_conditions)
    cond <- factor(new_conditions$conditions) %>% relevel(control)
    print(cond)
    design <- model.matrix(~ cond)
    fit <- lmFit(gsva_kegg[,new_conditions$new_samples], design1)
    fit<- eBayes(fit, trend=T)
    
    #colnames(exp_data) <- new_conditions$new_samples
    ## need to
  } else { #plot and analyse ALL the data (no restrictions)
    cond <- factor(conditions) %>% relevel(control_cond) # DMSO is the control
    design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
    colnames(design)[1] <- c(control_cond) 
    colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                               nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
    fit <- lmFit(gsva_kegg, design)
    fit <- eBayes(fit, trend=T)
    
    colnames(exp_data) <- new_samples
    new_conditions <- data.frame(new_samples, conditions)
  }
  
  list(exp_data = exp_data, fit = fit, design = design, gsva_kegg = gsva_kegg, conditions = conditions,
       new_conditions = new_conditions)
})
