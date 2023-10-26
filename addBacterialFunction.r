# function to add functional groups to bacterial taxonomy
# to create function

addBacterialFunction <- function(tax, tax_fun_ref = NULL){
  
  if (is.null(tax_fun_ref)) {
    cat("tax_fun_ref argument is empty. Use 'createTaxFunction()' to create the required dataset.")
    return()
  }
  
	# save rownames for later
	tax <- as(tax, "matrix")
	rowname <- rownames(tax)
	colnames(tax) <- tolower(colnames(tax))
	tax <- as.data.frame(tax)
	tax <- data.frame(lapply(tax, function(x) tolower(x)), stringsAsFactors = FALSE)
	
  # Classifications from literature search (multiple taxon levels)
  pathway_names <- colnames(tax_fun_ref)[3:ncol(tax_fun_ref)]
  tax[, pathway_names] <- "other"

  # taxon assignments
  for (i in 1:length(pathway_names)) {
    pathway <- pathway_names[i]
    levels <- c("phylum", "class", "order", "family", "genus")
    for (j in 1:length(levels)) {
      taxon_level <- levels[j]
      taxa_with_pathway <- tax_fun_ref %>% # get subset of taxa in group *at the correct level*
        dplyr::filter(tax_fun_ref[,!!pathway]==1) %>% 
        dplyr::filter(taxonomic.level == taxon_level) 
      if(nrow(taxa_with_pathway)==0) next()
      #print(taxon_level) #testing
      taxa_with_pathway <- taxa_with_pathway$taxon
      tax <- tax %>% # assign pathway as present for those taxa, and leave already-classified taxa alone
        dplyr::mutate(!!pathway := dplyr::case_when(tax[,!!taxon_level] %in% taxa_with_pathway ~ !!pathway,
                                      tax[,!!pathway] == !!pathway ~ !!pathway,
                                      TRUE ~ "other"))
    }
  }
  rownames(tax) <- rowname
  return(as.matrix(tax))
}
