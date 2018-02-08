library(jsonlite)
library(reshape2)
#library(igraph)

##### load data

# load protein - glycan - motif correspondence
#glycan_motif2 = fromJSON('data/glycan_motifs/match_full.json')
glycan_motif_protein = fromJSON('data/annotation/glycan_uniprotkbID_match_with_motif_hit.json')
glycan_protein = melt(lapply(glycan_motif_protein,function(x) x[[2]]))
colnames(glycan_protein) = c('uniprotswissprot','glytoucanID')
glycan_motif = melt(do.call(cbind,lapply(glycan_motif_protein,function(x) x[[1]])))
colnames(glycan_motif) = c('motif','glytoucanID','occurance')
#protein_glycan = fromJSON('data/glycan_motifs/uniCarbKB_Homo_glycan_dic.json')
#protein_glycan = fromJSON('data/annotation/Human_glycan_uniprotkbID_match.json')
#glycans = fromJSON('data/glycan_motifs/glycan_within_all_lectins.json')

# load functional annotation
GO = read.csv('data/annotation/annotated_genes.GO.txt')
#adj_GO = as_adjacency_matrix( graph_from_edgelist( as.matrix(GO[,2:3]) , directed=F))
#adj_GO = as.matrix(adj_GO[grepl('^GO:',rownames(adj_GO)),!grepl('^GO:',colnames(adj_GO))])
domains = read.csv('data/annotation/annotated_genes.interpro.txt')
domains[domains=='']=NA
domains=na.omit(domains)
#adj_dom = as_adjacency_matrix( graph_from_edgelist( as.matrix(domains[,2:3]) , directed=F))
#adj_dom = as.matrix(adj_dom[grepl('^IPR',rownames(adj_dom)),!grepl('^IPR',colnames(adj_dom))])

##### run stats

# regression model
library(doParallel)
library(foreach)
cl<-makeCluster(spec = 20)
registerDoParallel(cl = cl)
#for( m in unique( glycan_motif$motif))){
jnk=foreach( m=unique(glycan_motif$motif) , .combine=c) %dopar% {
	out=list()
	gmp_GO = merge( merge( glycan_protein , glycan_motif[glycan_motif$motif==m,] ) , GO[,2:3] ) 
	for( go in levels( gmp_GO$go_id)){
		data = data.frame( go_i = as.numeric(gmp_GO$go_id==go) , gmp_GO)
		# check assumtions: 
		tab = table( data$go_i , data$occurance )
		eo = apply(tab,1,function(x) apply(tab,2,function(y) (sum(x)*sum(y))/sum(tab)))
		if( any(eo<10) ){next} # cochran(1952,1954)
#		if( any(eo<1) | sum(eo<5)<(.2*prod(dim(tab))) ){next} # Yates, Moore & McCabe 1999 (tables larger and 2x2)
		prop = sum(data$go_i)/nrow(data)
		data$weights = ifelse( data$go_i==1 , 1 , prop)
		out[[paste0(go,'_motif',m)]] = glm( go_i ~ occurance , data=data,weights=data$weights,family='binomial')
#		out[[paste0(go,'_motif',m)]] = glmer( go_i ~ occurance + (1|uniprotswissprot) , data=data,weights=data$weights , family='binomial')
	}
	save(out,file=paste0('associations/',m,'.out.rda'))
	NULL
}
stopCluster(cl)

#hist( tmp<-unlist(lapply(out,function(x) coef(summary(x))[2,4] )))
