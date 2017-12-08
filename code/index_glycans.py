## index the glycans into a string searchable datastructure
# input: dictionary: {glycan ID:IUPAC strings}
# output: suffix tree

def index_glycans(G):
	IDX = SuffixTree()
	DFSlen = {}
	for k in G.keys():
		g = G[k]
		# append glycan index to the end of each glycan to make them terminally uniquely identifiable
		gk = g+":"+k
		# split glycan into every unique string by DFS
		gkL = glycanDFS(gk)
		DFSlen[k] = len(gkl)
		# add strings to suffix tree 
		for(s in gkL):
			IDX.add(s)
	return(IDX,DFSlen)
