## Example:

# Example data created from:
#example_data <- Allen_downsampled_logCPM_dataset %>%
#  dplyr::select(gene_symbol, L1)
#save(example_data, file = '/Users/ethankim/Desktop/example_data.Rdata')

#### Script for analysis

# Source functions
source(file='/Users/ethankim/Desktop/GO_enrichment_functions.R') # change for your own file dir
# Load data
load(file='/Users/ethankim/Desktop/example_data.Rdata') # change for your own dir

# Create background list of genes
gene_background_list <- create_geneBackground(example_data, 'gene_symbol')

# Create ranked list of genes
ranked_gene_list <- generate_AUC_genelist(example_data, 'L1', 'gene_symbol')

# Create geneSet tmod object
# This DOES take a long time - pls be patient :)
gene_sets_tmod_GO <- create_geneSets_GO(gene_background_list, 20, 100) # numbers are up to you!

# Create table of AUC values + other metadata
AUC_table <- create_AUC_table(ranked_gene_list, gene_sets_tmod_GO)

# Show table preview
head(AUC_table)

# Metadata:
# MainTitle: GO group name
# N1: number of genes within group
# AUC: AUC value
# P.value: P value associated with AUC
# adj.P.Val: P value adjusted by FDR
# ID: GO group ID
# aspect: GO aspects (CC: cellular component; MF: molecular function; BP: biological process)
# otherNames: other names for group
# rank: rank of AUC value, from lowest to highest value
