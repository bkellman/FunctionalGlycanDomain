library(broom)
# setwd('../Google Drive/00_professional/projects/FunctionalGlycanDomain/')

r = read.csv('data/glycan_motifs/lectinsVSmotifs.txt')

motifs = grep('X.G00',colnames(r),value = T)
taxa = c(' X.Genus.','X.Species.')

##
library(gplots)
heatmap.2( na.omit(data.matrix(r[,motifs])) , trace='none', labRow = r$X.GenomeComp.)
heatmap.2( na.omit(data.matrix(r[order(r$X.GenomeComp.),motifs])) , trace='none', labRow = r$X.GenomeComp.,Rowv = F)

## init - genometype variables
r$X.GenomeComp.[r$X.GenomeComp.=='ssrna-rt'] = NA
r$ss_ds = factor( ifelse( grepl('^ds',r$X.GenomeComp.) , 'ds', ifelse(grepl('^ss',r$X.GenomeComp.), 'ss', NA)))
r$rna_dna = factor( ifelse( grepl('rna',r$X.GenomeComp.) , 'rna', ifelse(grepl('dna',r$X.GenomeComp.), 'dna', NA)))
r$strand = factor( ifelse( grepl('rna\\(\\-\\)',r$X.GenomeComp.) , 'neg', ifelse(grepl('rna\\(\\+\\)',r$X.GenomeComp.), 'pos', NA)))

type = c('ss_ds','rna_dna','strand')

## associations
l=list()
for(t in type){
  for(m in motifs){
    if( length(table(r[[m]]))<2 ){next}
    tab=table(r[[t]],r[[m]],useNA = 'no')
    expected = chisq.test(tab)$expected
    attr(tab,'fisher_assumptions')= sum(tab)<1000 & all(expected>1) & (sum(expected<5)/dim(expected))>.2
    attr(tab,'chisqr_assumptions')= all(expected>1) & (sum(expected>5)/dim(expected))>.8
    l[[paste0(t,'_',m)]]=tab
  }
}

fisher=lapply(names(l),function(n){
  li=l[[n]]
  if( attr(li,'fisher_assumptions')){
    print(n)
    print(paste('Pr(fisher):',fisher.test(li)$p.value))
    print(li)
  }
})

chisqr=lapply(names(l),function(n){
  li=l[[n]]
  if( attr(li,'chisqr_assumptions')){
    print(n)
    print(paste('Pr(chisqr):',chisq.test(li)$p.value))
    print(li)
  }
})


