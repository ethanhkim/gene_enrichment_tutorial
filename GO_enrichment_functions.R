# Load required packages

require(tmod)
require(dplyr)
require(magrittr)
require(org.Hs.eg.db)
require(tibble)
require(annotate)
require(AnnotationDbi)
require(GO.db)
library(docstring) # for python-esque docstring - optional

#### Functions


create_geneBackground <- function(source_data, column_of_gene_symbols) {
  
  #' @description
  #' Creates a named list of the full list of genes to use in your AUC analysis
  
  #' * source_data: dataset you're working with
  #' * column_of_gene_symbols: string of column of the full list of genes you want to 
  #' * use as the background for the AUC analysis
  
  # Pull the column of gene symbols
  geneBackground <- source_data %>%
    # Select gene_symbol column
    dplyr::select(!!column_of_gene_symbols) %>%
    # Pull first column
    pull()
  # Set names
  names(geneBackground) <- geneBackground
  # Sort alphabetically, and only keep duplicates
  geneBackground <- sort(unique(names(as.list(geneBackground))))
  return(geneBackground)
}

# Create geneSetsTable for analysis for Gene Ontology annotations
create_geneSets_GO <- function(geneBackground_list, min_GO_group_size, max_GO_group_size) {
  
  #' @description
  #' Creates a geneSet tmod object used for AUC analysis containing GO groups
  
  #' * geneBackground_list: list of genes created from create_geneBackground
  #' * min_GO_group_size: minimum number of genes to be contained in GO group
  #' * max_GO_group_size: maximum number of gnees to be contained in GO group
  
  # Create GO object
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  # Get gene_symbols in GO_object
  symbols_in_GO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
    
  #build GO sets for tmod -slow
  tmod_names <- data.frame()
  modules2genes <- list()
  go_group_name <- names(go_object)[1]
  showMethods(Term)
    
  go_count <- length(go_object)
  count <- 1
  for(go_group_name in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", go_count))
    count <- count + 1
    # Select by go_group_name
    go_group <- go_object[go_group_name]
    gene_IDs <- unique(unlist(go_group, use.names=F))  #discard evidence codes
    gene_symbols <- unique(getSYMBOL(gene_IDs, data='org.Hs.eg'))
    #get size after intersecting with our full gene set
    genesymbols <- intersect(gene_symbols, geneBackground_list) 
    if (!(length(gene_symbols) >= min_GO_group_size & length(gene_symbols) <= max_GO_group_size)) next();
      
    modules2genes[go_group_name] <- list(gene_symbols)
    tmod_names <- rbind(tmod_names, data.frame(ID=go_group_name, Title = Term(go_group_name)))
  }
  # Create tmod object from GO sets
  geneSetsGO <- makeTmod(modules = tmod_names, modules2genes = modules2genes)
  # Return geneSetsGO tmod object
  return(geneSetsGO)
}


generate_AUC_genelist <- function(source_data, column_to_order, column_of_gene_symbols) {
  
  #' @description
  #' Creates an ordered list of genes arranged by descending order of expression
  
  #' * source_data: dataset you're working with
  #' * column_of_gene_symbols: column of the full list of genes 
  #' * column_to_order: column in data containing gene expression values to rank
  
  # Select coluns specified
  selected_data <- source_data %>%
    dplyr::select(!!column_of_gene_symbols, !!column_to_order)
  # Sort by expression is descending order
  sorted_data <- selected_data[order(selected_data[, column_to_order], decreasing = TRUE),]
  # Pull gene symbol column, which should be sorted by expression
  genelist <- sorted_data %>% pull(!!column_of_gene_symbols) %>% unique()
  # Return genelist
  return(genelist)
}



create_AUC_table <- function(AUC_genelist, geneSet) {
  
  #' @description
  #' Creates table of AUC values, and associated statistics from tmod
  
  #' * AUC_genelist: genelist from generate_AUC_genelist function
  #' * geneSet: geneSetsGO from create_geneSets_GO function
  
  # Creates a top and bottom 20 list of ontology definitions using AUC
  AUC_table <- as_tibble(tmodUtest(c(AUC_genelist), mset=geneSet, qval = 1.01, filter = T))
  AUC_table %<>% 
    rowwise() %>% 
    #tmod runs one-sided tests - remove if performing 1-sided test
    mutate(P.Value = P.Value * 2) %>% 
    ungroup() %>% 
    # Create adjusted p-value column, adjusted by FDR
    mutate(adj.P.Val=p.adjust(P.Value, method="fdr")) 
  AUC_table %<>% rowwise() %>% mutate(aspect = Ontology(ID))
  
  # Collapse genesets that have the exact same set of genes
  AUC_table %<>% rowwise() %>% 
    mutate(genes = paste(sort(unlist(geneSet$MODULES2GENES[ID])), collapse = " "))
  AUC_table %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = dplyr::first(Title),  ID=paste(ID, collapse=","), AUC = dplyr::first(AUC), 
              P.Value= dplyr::first(P.Value), aspect =dplyr::first(aspect), 
              otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  
  AUC_table %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="fdr"))
  AUC_table %<>% dplyr::select(-genes)
  AUC_table %<>% arrange(P.Value)
  # Create rank of AUC values
  AUC_table$rank <- 1:nrow(AUC_table)
  # Select columns
  AUC_table %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)
  # Return table
  return(AUC_table)
}

