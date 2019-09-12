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
match_pathway <- function(df, annot_type, core_pep_kegg){
  
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