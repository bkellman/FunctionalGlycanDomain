from biomart import BiomartServer

def GOquery(query,organism,gene_id):
	# connect to server
	server = BiomartServer( "http://www.biomart.org/biomart" )
	# set verbose to True to get some messages
	server.verbose = True

	db = server.datasets[organism+'_gene_ensembl']

	response = db.search({
	'filters': { gene_id: query },
	'attributes': [ gene_id, 'with_go' ]
	})
	return response
