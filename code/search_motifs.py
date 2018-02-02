## search motifs
# input:
#    dictionary of glycan motifs in IUPAC form
#    Index of glycans
# output:
#    Dictionary with glycans ID as keys matched to a list of all glycan motifs they contain


def search_motifs(M,IDX):
    Map={}
    for k in M.keys():
        m = M[k]
        # slit into substrings
        S = glycanDFS(m)
        # find each string in the index mapping to string s of motif m
        match = [IDX.find(s) for s in S]
        # Identify all glycans contain EVERY DFS substring s of motif m
        match_allGlycans_bySubstring = [mi.split(':') for mi in match]
        match_uniqueGlycans = Set(flatten( match_allGlycans_bySubstring ))
            # what proportion of the substrings for motif m are discover glycan g?
        prop = [ {g:sum(inL(g,match_allGlycans_bySubstring))/len(S)} for g in match_uniqueGlycans]
        Map[k] = []
        for g in prop.keys():
            if prop[g]==1: # if 100% of the motif substrings are found within glycan g, include glycan g as mapping to motif k
                Map[k].append(g)
    return(Map)


def flatten(l):
    [item for sublist in l for item in sublist]

def inL(x,L):
    out=[]
    for l in L:
        if x in l:
            out.append(1)
        else:
            out.append(0)
    return(out)
