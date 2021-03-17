#' Generate list of relevant probes for gain, loss, and change in methylation 
#' status 
#'
#' This function returns relevant probe lists for each node that exhibits 
#' either a gain or a loss in methylation status. It also returns three sets 
#' of branch labels: number of methylation changes, number of methylation 
#' gains, and number of methylation losses. These labels can be given to the 
#' show_tree function as an edge labels parameter to visualize this data.
#' 
#' @param icceTree icceTree data structure
#' @param probe_node_matrix probe-node matrix storing the methylation status 
#' across all probes and nodes on a tree
#'
#' @return a structure containing a list of probes with differentially 
#' methylated probes for each node (gain and loss stored separately), a list 
#' of branch labels for the number of total changes in methylation, a list of 
#' branch labels for the number of gains in methylation, and a list of branch 
#' labels for the number of losses in methylation
#' @examples
#' probe_information <- generate_branch_changes(icceTree)
#' relevant_probes <- probe_information$relevant_probes
#' total_changes_label <- probe_information$total_changes_label
#' meth_gains_label <- probe_information$meth_gains_label
#' meth_losses_label <- probe_information$meth_losses_label
#' @export
generate_branch_changes <- function(icceTree) {
	tree <- icceTree$tree
	probe_node_matrix <- icceTree$probe_node_matrix
	branch_changes = data.frame(total_changes=integer(), meth_gains=integer(), meth_losses=integer())

	relevant_probes = c()

	for (i in 1:nnodes(tree)) {

		row = as.matrix(tree$edge[i,])
		i_ancestor = toString(row[1])
		i_descendant = toString(row[2])

		probe_node_matrixcol = rownames(probe_node_matrix)

		branch_changes[i, 'total_changes'] = length(probe_node_matrixcol[probe_node_matrix[[i_ancestor]] != probe_node_matrix[[i_descendant]]])

		meth_gains_probes = probe_node_matrixcol[probe_node_matrix[[i_ancestor]] != probe_node_matrix[[i_descendant]] & probe_node_matrix[[i_descendant]] == 1 & probe_node_matrix[[i_ancestor]] == 0]
		branch_changes[i, 'meth_gains'] = length(meth_gains_probes)

		meth_losses_probes = probe_node_matrixcol[probe_node_matrix[[i_ancestor]] != probe_node_matrix[[i_descendant]] & probe_node_matrix[[i_descendant]] == 0 & probe_node_matrix[[i_ancestor]] == 1]
		branch_changes[i, 'meth_losses'] = length(meth_losses_probes)

		relevant_probes[[paste(toString(i), "_gain", sep = "")]] = meth_gains_probes
		relevant_probes[[paste(toString(i), "_loss", sep = "")]] = meth_losses_probes
	}

	# Define the branch labels based on the total number of changes
	total_changes_label = as.numeric(branch_changes$total_changes)

	# Define branch labels based on the number of methyl gains
	meth_gains_label = as.numeric(branch_changes$meth_gains)

	# Define branch labels based on the number of methyl losses
	meth_losses_label = as.numeric(branch_changes$meth_losses)

	return(list('relevant_probes' = relevant_probes, 'total_changes_label' = total_changes_label, 'meth_gains_label' = meth_gains_label, 'meth_losses_label' = meth_losses_label))
}

#' Find relevant genes for gain and loss of methylation 
#'
#' This function returns a list of relevant gene for each node that exhibit a 
#' gain or loss in methylation status. 
#' 
#' @param icceTree icceTree data structure
#' @param relevant_probes list of differentially methylated probes for each 
#' node
#' @param reference reference file
#'
#' @return list of differentially methylated genes for each node
#' @examples
#' relevant_genes <- get_relevant_genes(icceTree, relevant_probes, reference)
#' any(grepl(gene, bcells))... can be used, create function for this
#' @export
get_relevant_genes <- function(icceTree, relevant_probes, reference) {
	tree <- icceTree$tree

	relevant_genes = c()

	edges = paste(tree$edge[,1], tree$edge[,2], sep="-")
	gain_name = paste(edges, "_gain", sep="")
	loss_name = paste(edges, "_loss", sep="")
	edge_names = append(gain_name, loss_name)

	for (i in 1:nnodes(tree)) {
		i_gain = paste(toString(i), "_gain", sep = "")
		i_loss = paste(toString(i), "_loss", sep = "")

		# Retrieve genes from specified probes
		relevant_genes[[i_gain]] = reference[rownames(reference) %in% relevant_probes[[i_gain]], "gene"]
		relevant_genes[[i_loss]] = reference[rownames(reference) %in% relevant_probes[[i_loss]], "gene"]

		# Remove NA genes
		relevant_genes[[i_gain]] = relevant_genes[[i_gain]][!is.na(relevant_genes[[i_gain]])]
		relevant_genes[[i_loss]] = relevant_genes[[i_loss]][!is.na(relevant_genes[[i_loss]])]

		i_descendant = as.numeric(tree$edge[i, 1])
		i_ancestor = as.numeric(tree$edge[i, 2])

		names(relevant_genes)[2 * i - 1] = edge_names[grepl(i_descendant, edge_names) & grepl(i_ancestor, edge_names) & grepl("gain", edge_names)]

		names(relevant_genes)[2 * i] = edge_names[grepl(i_descendant, edge_names) & grepl(i_ancestor, edge_names) & grepl("loss", edge_names)]

	}

	return(relevant_genes)
}

#' Get probes along a given gene 
#'
#' This function returns a list of probes along a particular gene using the 
#' given reference file.
#' 
#' @param reference reference file
#' @param gene gene name
#'
#' @return list of probes along a particular gene
#' @examples
#' gene_probes <- get_probes_of_gene(reference, gene)
#' @export
get_probes_of_gene <- function(reference, gene) {
	gene_probes = rownames(reference[reference$gene == gene & !is.na(reference$gene), ])
	return(gene_probes)
}

#' Get proportion of probes methylated, unmethylated, ambiguous at each node 
#' across a tree
#'
#' This function returns the proportion of each probes that are methylated, 
#' unmethylated, ambiguous at each node across a tree.
#' 
#' @param icceTree icceTree data structure
#' @param probe_node_matrix matrix storing the methylation status across all
#' probes and nodes on a tree
#' @param reference reference file
#' @param gene gene name
#'
#' @return a structure composed of two lists that contain methylation state 
#' proportions for internal nodes and methylation state proportions for leaves
#' @examples
#' gene_information <- track_gene(icceTree, reference, 'CD81')
#' thermo_prop_internal_nodes = gene_information$thermo_prop_internal_nodes
#' thermo_prop_leaves = gene_information$thermo_prop_leaves
#' @export
track_gene <- function(icceTree, reference, gene) {
	tree <- icceTree$tree
	probe_node_matrix <- icceTree$probe_node_matrix

	gene_probes <- get_probes_of_gene(reference, gene)

	probe_node_matrix_states = probe_node_matrix[rownames(probe_node_matrix) %in% gene_probes,]

	probe_node_matrix_states = data.frame(t(probe_node_matrix_states))

	probe_node_matrix_states[probe_node_matrix_states == "1.0"] = 1
	probe_node_matrix_states[probe_node_matrix_states == "0.0"] = 0

	nprobes = ncol(probe_node_matrix_states)

	thermo_prop = c()

	thermo_prop$"1" = rowSums(probe_node_matrix_states[,1:nprobes] == 1) / nprobes

	thermo_prop$"0" = rowSums(probe_node_matrix_states[,1:nprobes] == 0) / nprobes

	thermo_prop$"0.5" = rowSums(probe_node_matrix_states[,1:nprobes] == 0.5) / nprobes

	thermo_prop = data.frame(thermo_prop, check.names=FALSE, row.names = NULL)

	rownames(thermo_prop) = rownames(probe_node_matrix_states)

	# Use leaves(tree) function

	thermo_prop_internal_nodes <- thermo_prop[1:5,]

	thermo_prop_leaves <- thermo_prop[6:11,]

	return(list('thermo_prop_internal_nodes' = thermo_prop_internal_nodes, 'thermo_prop_leaves' = thermo_prop_leaves))
}
