###########
## load glycans
#######

experiment = data.frame( factor1 = rep(letters[1:2],each=2) ,factor2 = rep(letters[1:2],2) )

abundance_g = matrix( runif(20) , nrow=5,ncol=4)
rownames( abundance_g) = paste('glycan',1:5)
colnames( abundance_g) = paste('profile',1:4)

###########
## extranct glycan motifs and motif relations
#######

# for each glycan
	extract glycan motifs

# motif x motif relation (parent in rows, child in columns)
relation = matrix( 0 , 10,10)
rownames( relation) = paste('motif',1:10)
colnames( relation) = paste('glycan',1:5)

# for each motif
	# for each motif
		# is motif1 and member of motif2?
			# if so, relations[motif1,motif2] = 1

# reduce redundant edges

###########
## motif quantification
#######

# initialize motif abunance
abundance_m = matrix( 0 , nrow=10,ncol=4)
rownames( abundance_m) = paste('motif',1:10)
colnames( abundance_m) = paste('profile',1:4)

## calculate motif abunace
for each glycan
	for each motif
		if motif is in a glycan,
			abundance_m[motif,] = abundance_m[motif,] + abundance_g[glycan,] 

###########
#### THE STATS 
############

#### background agnostic
for( f in 1:ncol(experiment)){
	f_i = experiment[,f]
	for( m in motifs)
		glm ( f_i ~ m )  
}

#### background consideration (LRT)
parent_of <- function(m1,m2,relation){
	return all nonzero (relations[,m2])
}

for( f in 1:ncol(experiment)){
	f_i = experiment[,f]
	for( m in motifs)
		null_model = glm ( f_i ~ parents_of(m) )  
		alt_model = glm( f_i ~ parents_of(m) + m) 
		anova( null_model , alt_model , test='LRT') 
}
