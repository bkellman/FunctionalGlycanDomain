# load glycans
G = loadGlycans('...')

M = loadMotifs('...')

# assemble index
idx = index_glycans(G)

# find motifs in glycans
Dg2m = search_motifs(G,M)
