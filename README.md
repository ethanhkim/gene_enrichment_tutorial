# Welcome to the GO enrichment analysis tutorial!

This tutorial repo is meant to provide a quick tutorial on how to derive the enrichment of Gene Ontology (GO) terms in a given list of genes and their expression.

### What's included in this handy dandy tutorial?

* `GO_enrichment_functions.R`: the functions necessary for examining GO enrichment
* `GO_enrichment_script.R`: an example script
* `example_data.Rdata`: an example dataset of genes and their respective expression

### Included functions:

* `create_geneBackground`: creates a list of   
* `create_geneSets_GO`: creates a tmod object made up of GO objects that include gene symbols in HGNC
* `generate_AUC_genelist`: ranks genes by their expression
* `create_AUC_table`: uses tmod object created by `create_geneSets_GO` and p

### Instructions:

1. Pull the full list of genes to create the background using `create_geneBackground`
2. Rank the dataset by using your genes' expression with `generate_AUC_genelist`
3. Create GO sets of your full list of genes with `create_geneSets_GO` 
4. Examine the enrichment of GO terms from your own gene list with `create_AUC_table`

An example script (`GO_enrichment_script.R`) can be found to reference!