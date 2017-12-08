## Depth First Search Linearization of an IUPAC glycan string
# Input: glycan tree (in linear IUPAC form)
# Output: List of all strings discoverable by depth first search

import pygraph

def search_motifs(IUPAC):
	G = IUPAC2graph(IUPAC)
	paths = find_all_paths(G,G.root,G.leafs)
	return(paths)

def find_all_paths(graph, start, end, path=[]):
	path = path + [start]
	if start == end:
		return [path]
	if not graph.has_key(start):
		return []
	paths = []
	for node in graph[start]:
		if node not in path:
			newpaths = find_all_paths(graph, node, end, path)
			for newpath in newpaths:
				paths.append(newpath)
	return paths
