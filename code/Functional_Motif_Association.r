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

run_go=T
run_domains=T

if(run_go){
	#### run GO
	# regression model
	library(doParallel)
	library(foreach)
	cl<-makeCluster(spec = 20)
	registerDoParallel(cl = cl)
	#for( m in unique( glycan_motif$motif))){
	jnk=foreach( m=unique(glycan_motif$motif) ,.errorhandling='pass' ) %dopar% {
		out=list()
		if( file.exists(paste0('associations/',m,'.out.rda')) ) { NULL }
		gmp_GO = merge( merge( glycan_protein , glycan_motif[glycan_motif$motif==m,] ) , GO[,2:3] ) 
		for( go in levels( gmp_GO$go_id)){
			data = data.frame( go_i = as.numeric(gmp_GO$go_id==go) , gmp_GO)
			# check assumtions: 
			tab = table( data$go_i , data$occurance )
			if(min(dim(tab))<2){next}
			eo = apply(tab,1,function(x){ apply(tab,2,function(y){
				x=as.numeric(x); y=as.numeric(y); tb = as.numeric(tab)
				(sum(x)*sum(y))/sum(tb)
			})})			
			if( any(eo<10) ){next} # cochran(1952,1954)
	#		if( any(eo<1) | sum(eo<5)<(.2*prod(dim(tab))) ){next} # Yates, Moore & McCabe 1999 (tables larger and 2x2)
			prop = sum(data$go_i)/nrow(data)
			data$weights = ifelse( data$go_i==1 , 1 , prop)
			out[[paste0(go,'_motif',m)]] = glm( go_i ~ occurance , data=data,weights=data$weights,family='binomial')
	#		out[[paste0(go,'_motif',m)]] = glmer( go_i ~ occurance + (1|uniprotswissprot) , data=data,weights=data$weights , family='binomial')
		}
		save(out,file=paste0('associations/',m,'.out.rda'))
		gc(reset=T)
		NULL
	}
	stopCluster(cl)

	entry_count = list()
	out=do.call(rbind,tmp<-lapply( system('ls associations/*.out.rda',inter=T) , function(file){
		print(file)
		i = strsplit(file,'\\.')[[1]][1]
		#load(paste0('associations/',file))
		load(file)
		if(length(out)>0){
			outi=as.data.frame(do.call(rbind,lapply(out,function(o){
				coef(summary(o))[2,]
			})))
			outi$comparison = names(out)
		}else{outi=NULL}
		entry_count[[as.character(i)]] = length(out)
		outi
	}))
	save(out,file='associations/go_out.rda')
	system('rm associations/*.out.rda')
}

if(run_domains){
	#### run domains
	# regression model
	library(doParallel)
	library(foreach)
	cl<-makeCluster(spec = 20)
	registerDoParallel(cl = cl)
	#for( m in unique( glycan_motif$motif))){
	jnk=foreach( m=unique(glycan_motif$motif) ,.errorhandling='pass' ) %dopar% {
		out=list()
		if( file.exists(paste0('associations/',m,'.out.rda')) ) { NULL }
		gmp_domains = merge( merge( glycan_protein , glycan_motif[glycan_motif$motif==m,] ) , domains[,3:4] ) 
		for( dom in levels( gmp_domains$interpro)){
			data = data.frame( dom_i = as.numeric(gmp_domains$interpro==dom) , gmp_domains)
			# check assumtions: 
			tab = table( data$dom_i , data$occurance )
			if(min(dim(tab))<2){next}
			eo = apply(tab,1,function(x){ apply(tab,2,function(y){
				x=as.numeric(x); y=as.numeric(y); tb = as.numeric(tab)
				(sum(x)*sum(y))/sum(tb)
			})})
			if( any(eo<10) ){next} # cochran(1952,1954)
	#		if( any(eo<1) | sum(eo<5)<(.2*prod(dim(tab))) ){next} # Yates, Moore & McCabe 1999 (tables larger and 2x2)
			prop = sum(data$dom_i)/nrow(data)
			data$weights = ifelse( data$dom_i==1 , 1 , prop)
			out[[paste0(dom,'_motif',m)]] = glm( dom_i ~ occurance , data=data,weights=data$weights,family='binomial')
	#		out[[paste0(go,'_motif',m)]] = glmer( go_i ~ occurance + (1|uniprotswissprot) , data=data,weights=data$weights , family='binomial')
		}
		save(out,file=paste0('associations/',m,'.out.rda'))
		gc(reset=T)
		NULL
	}
	stopCluster(cl)

	entry_count = list()
	out=do.call(rbind,tmp<-lapply( system('ls associations/*.out.rda',inter=T) , function(file){
		print(file)
		i = strsplit(file,'\\.')[[1]][1]
		# load(paste0('associations/',file))
		load(filed)
		if(length(out)>0){
			outi=as.data.frame(do.call(rbind,lapply(out,function(o){
				coef(summary(o))[2,]
			})))
			outi$comparison = names(out)
		}else{outi=NULL}
		entry_count[[as.character(i)]] = length(out)
		outi
	}))
	save(out,file='associations/domain_out.rda')
	system('rm associations/*.out.rda')
}

# ## vis
# load('associations/out.rda')
# library(RamiGO)

# p=.05
# q=.1
# e=1.5
# out_1860_neg = unlist(lapply(strsplit(out$comparison[grepl('motif1860',out$comparison) & p.adjust(out[['Pr(>|z|)']],'fdr')<.1 & out$Estimate<(-e)],'_'),function(x) x[1]))
# out_1860_pos = unlist(lapply(strsplit(out$comparison[grepl('motif1860',out$comparison) & p.adjust(out[['Pr(>|z|)']],'fdr')<.1 & out$Estimate>e],'_'),function(x) x[1]))

# out_1840_neg = unlist(lapply(strsplit(out$comparison[grepl('motif1840',out$comparison) & p.adjust(out[['Pr(>|z|)']],'fdr')<.1 & out$Estimate<(-e)],'_'),function(x) x[1]))
# out_1840_pos = unlist(lapply(strsplit(out$comparison[grepl('motif1840',out$comparison) & p.adjust(out[['Pr(>|z|)']],'fdr')<.1 & out$Estimate>e],'_'),function(x) x[1]))

# goIDs <- unique( c(out_1840_neg,out_1840_pos) )
# tab = table(goIDs)
# color <- ifelse(goIDs %in% out_1840_neg,
# 	ifelse(goIDs %in% out_1840_pos,'green','yellow'),
# 	ifelse(goIDs %in% out_1840_pos,'blue',NA))
# pngRes <- getAmigoTree(goIDs=goIDs, color=color, filename="m1840_neg_pos", picType="png", saveResult=TRUE)