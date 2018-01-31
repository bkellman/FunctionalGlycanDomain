library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

query = read.csv('gene_names.txt')$uniprotKB_id

out <- getBM(attributes=c('uniprot_gn','hgnc_symbol','uniprotswissprot','name_1006','definition_1006'), 
	filters = 'uniprotswissprot', values =query, mart = ensembl)
	
write.csv(out,file='annotated_genes.txt')