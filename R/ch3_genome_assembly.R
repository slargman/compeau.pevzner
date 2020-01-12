#' Generate the k-mer composition of a string.
#' 
#' Given a string \code{text}, its \emph{k}-mer composition is the collection of all \emph{k}-mer substrings of \code{text} (including repeated \emph{k}-mers. \code{StringComposition generates the \emph{k}-mer composition of a given string, listed in lexicographic order.
#' 
#' @param text A string (i.e. a DNA sequence consisting only of the characters A, C, G, or T)
#' @param k An integer less than or equal to the number of characters in \code{text}
#' @return A character vector listing all \emph{k}-mer substrings of \code{text} in lexicographic order.
#' @examples
#' StringComposition("TATGGGGTGC", 3)
#' StringComposition("CAATCCAAC", 5)
StringComposition <- function(text, k){
	n <- nchar(text)
	kmer_composition <- substring(text, 1:(n - k + 1), k:n)
	kmer_composition <- sort(kmer_composition)
	return(kmer_composition)
}

#' Reconstruct a string from its genome path
#' 
#' \code{StringFromGenomePath} reconstructs a string from a genome path of consecutive \emph{k}-mers that are listed in the character vector \code{pattern}. These consecutive \emph{k}-mers overlap such that the last \emph{k} - 1 characters of one element are the same as the first \emph{k} - 1 characters of the next. Concatenating these \emph{k}-mers yields the string from which they originate.
#' 
#' @param pattern A character vector of length \code{n} containing \emph{k}-mers such that the last \emph{k - 1} symbols of \code{pattern[i]} are equal to the first \emph{k - 1} symbols of \code{pattern[i + 1]} for \code{i} ranging from 1 to \emph{n} - 1.
#' @return A string of \code{text} of length \emph{k} + \emph{n} - 1 such that the \emph{i}-th \emph{k}-mer in \code{text} is equal to \code{pattern[i]} for \emph{i} ranging from 1 to \emph{n}.
#' @examples
#' pattern <- c("ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT")
#' StringFromGenomePath(pattern)
StringFromGenomePath <- function(pattern){
	# can't iterate from 2 for length == 1
	if (length(pattern) == 1){
		return(pattern)
	}
	text <- pattern[1]
	n <- length(pattern)
	k <- nchar(pattern[1])
	for (i in 2:n) {
		text <- paste0(text, substring(pattern[i], k))
	}
	return(text)
}

#' Extract the first \emph{k} - 1 characters of a \emph{k}-mer
#'
#' \code{Prefix} extracts the first \emph{k} - 1 characters of the \emph{k}-mer \code{pattern}.
#'
#' \code{Prefix} is used in the implementation of \code{\link{OverlapGraph}}.
#'
#' @param pattern A string.
#' @return A string consisting of all but the last character of \code{pattern}.
#' @examples
#' Prefix("AAT")
Prefix <- function(pattern){
	if (any(nchar(pattern) == 1)) {
		stop("pattern must have more than 1 character")
	}
	k <- nchar(pattern)
	prefix <- substring(pattern, 1, k - 1)
	return(prefix)
}

#' Extract the last \emph{k} - 1 characters of a \emph{k}-mer
#'
#' \code{Suffix} extracts the last \emph{k} - 1 characters of the \emph{k}-mer \code{pattern}.
#'
#' \code{Suffix} is used in the implementation of \code{\link{OverlapGraph}}.
#'
#' @param pattern A string.
#' @return A string consisting of all but the first character of \code{pattern}.
#' @examples
#' Prefix("TAA")
Suffix <- function(pattern){
	if (any(nchar(pattern) == 1)) {
		stop("pattern must have more than 1 character")
	}
	k <- nchar(pattern)
	suffix <- substring(pattern, 2, k)
	return(suffix)
}

#' Construct the overlap graph of a collection of strings
#' 
#' \code{OverlapGraph} generates the overlap graph of the collection of \emph{k}-mers in \code{pattern}, which is a directed graph. The nodes of this graph consist of all the \emph{k}-mers in \code{pattern}. An edge exists from \code{node1} to \code{node2} if the last \emph{k} - 1 characters of \code{node1} are the same as the first \emph{k} - 1 characters of \code{node2}, i.e. the \code{\link{Prefix}} of \code{node1} is equal to the \code{\link{Suffix}} of \code{node2}.
#' 
#' @param pattern A character vector where each element has the same number of characters.
#' @return A list representing the overlap graph of the collection of \emph{k}-mers \code{pattern}. This list has two named elements. The element \code{graph\$nodes} contains a character vector listing the nodes of the graph. The element \code{graph\$edges} contains a two column character matrix where each row corresponds to an edge in the graph and the columns list the two nodes of an edge.
#' @examples
#' pattern <- c("ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT")
#' OverlapGraph(pattern)
OverlapGraph <- function(pattern){
	n <- length(pattern)
	node1 <- numeric(0)
	node2 <- numeric(0)
	prefixes <- Prefix(pattern)
	for (node in pattern) {
		suffix <- Suffix(node)
		overlaps <- pattern[suffix == prefixes]
		node1 <- c(node1, rep(node, length(overlaps)))
		node2 <- c(node2, overlaps)
	}
	nodes <- sort(pattern)
	edges <- matrix(c(node1, node2), ncol = 2)
	colnames(edges) <- c("node1", "node2")
	lexicographic <- order(node1, node2)
	edges <- edges[lexicographic, ]
	graph <- list(nodes = nodes, edges = edges)
	return(graph)
}

#' Print the edges of a graph for input into Rosalind
#' 
#' \code{PrintEdges} returns a more readable form of the edges of a graph in the format necessary for Rosalind.
#' 
#' @param graph A (directed) graph in list form as generated by e.g.\code{\link{OverlapGraph}} or \code{\link{PathGraph}}. The graph consists of two named elements, a character vector \code{nodes} containing the nodes and a two-column character matrix \code{edges} containing the edges. The first column of the matrix \code{graph\$edges} contains the first node of an edge and the second column contains the second, with each row thus representing an edge.
#' @return A character vector where each element lists all the edges emanating from a node.
#' @examples
#' pattern <- c("ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT")
#' graph <- OverlapGraph(pattern)
#' PrintEdges(graph)
PrintEdges <- function(graph){
	display <- unique(graph$edges[, "node1"])
	for (i in seq_along(display)) {
		edge_index <- which(graph$edges[, "node1"] == display[i])
		node2 <- paste0(graph$edges[edge_index, "node2"], collapse = ",")
		display[i] <- paste(display[i], "->", node2)
	}
	return(display)
}

#' Generate the path graph for a string
#' 
#' \code{PathGraph} generates the \code{k}-mer path graph for a string \code{text} and specified \code{k}. Given a string \code{text}, the path graph of \code{k}-mers for it is the path consisting of \code{nchar(text)} - \code{k} + 1 edges, where the \emph{i}-th edge of this path is labeled by the \emph{i}-th \code{k}-mer in \code{text} and the \emph{i}-th node of the path is labeled by the \emph{i}-th (\code{k} - 1)-mer in \code{text}.
#' 
#' @inheritParams StringComposition
#' @return A list representing the path graph of the collection of \emph{k}-mers \code{pattern}. This list has two named elements. The element \code{graph\$nodes} contains a character vector listing the nodes of the graph. The element \code{graph\$edges} contains a three column character matrix where each row corresponds to an edge in the graph and the columns list respectively the outgoing node, the edge label, and the incoming node of an edge.
#' @examples
#' text <- "TAATGCCATGGGATGTT"
#' k <- 3
#' PathGraph(text, k)
PathGraph <- function(text, k){
	n <- nchar(text)
	# edges labelled by k-mers
	edge_labels <- substring(text, 1:(n - k + 1), k:n)
	# nodes labelled by (k - 1)-mers
	nodes <- substring(text, 1:(n - k + 2), (k - 1):n)
	m <- length(nodes)
	edges <- matrix(c(nodes[1:(m - 1)], edge_labels, nodes[2:m]), ncol = 3)
	colnames(edges) <- c("node1", "label", "node2")
	lexicographic <- order(edges[, "node1"], edges[, "node2"])
	edges <- edges[lexicographic, ]
	graph <- list(nodes = nodes, edges = edges)
	return(graph)
}

#' Construct the De Bruijn graph for a string or k-mer collection
#' 
#' DeBruijnGraph generates the De Bruijn graph for either a string or a collection of strings of equal length.
#' 
#' In order to generate the De Bruijn graph for a string \code{text}, an integer \code{k} must be provided. In this case the De Bruijn graph is constructed by generating the path graph of \code{k}-mers for the string using \code{\link{PathGraph} and gluing together identical nodes (i.e. \code{k}-mers).
#' 
#' If \code{k} is not provided then the DeBruijnGraph will interpret \code{text} as a collection of \code{k}-mers (thus all the elements of \code{text} should have the same length) and generate the De Bruijn graph. In this case the nodes of the graph are all the unique (\code{k} - 1)-mers occurring as a prefix or suffix of the element of \code{text} (see \code{\link{Prefix}} and \code{\link{Suffix}}. For each \code{k}-mer in \code{text}, we connect its prefix node to its suffix node by a directed edge, and use the \code{k}-mer as the edge label.
#' 
#' @param text A string or a character vector where each element has the same number of characters.
#' @inheritParams StringComposition
#' @return A list representing the DeBruijn graph of the single string or collection of \emph{k}-mers \code{text}. This list has two named elements. The element \code{graph\$nodes} contains a character vector listing the nodes of the graph. The element \code{graph\$edges} contains a three column character matrix where each row corresponds to an edge in the graph and the columns list respectively the outgoing node, the edge label, and the incoming node of an edge.
#' @examples
#' k <- 4
#' text <- "AAGATTCTCTAC"
#' DeBruijnGraph(text, k)
#' pattern <- c("GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG")
#' DeBruijnGraph(pattern)
DeBruijnGraph <- function(text, k = NULL){
	# check if constructing from k-mers
	if (is.null(k)) {
		prefix <- Prefix(text)
		suffix <- Suffix(text) 
		nodes <- sort(unique(c(prefix, suffix)))
		edges <- matrix(c(prefix, text, suffix), ncol = 3)
		colnames(edges) <- c("node1", "label", "node2")
		lexicographic <- order(prefix, suffix)
		edges <- edges[lexicographic, , drop = FALSE]
		graph <- list(nodes = nodes, edges = edges)
		return(graph)
	} else {
		graph <- PathGraph(text, k)
		graph$nodes <- unique(graph$nodes)
		return(graph)
	}
}

#' Convert an adjacency list into a graph
#' 
#' \code{AdjacencyListToGraph} converts a character vector containing an adjacency list to a graph. The adjacency list is the format used to specify a graph on Rosalind and the output graph generated by \code{AdjacencyListToGraph} matches the format used in \code{\link{OverlapGraph}} and \code{\link{PathGraph}}.
#' 
#' @param adj_list A character vector where each element gives the nodes adjacent to a particular node. The format for each entry is the name of the node followed by " -> " then a list of the adjacent nodes separated by commas.
#' @return A list representing the corresponding graph for the adjacency list \code{adj_list}. This list has two named elements. The element \code{graph\$nodes} contains a character vector listing the nodes of the graph. The element \code{graph\$edges} contains a two column character matrix where each row corresponds to an edge in the graph and the columns list the two nodes of an edge.
#' @examples
#' adj_list <- c("0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6") 
#' AdjacencyListToGraph(adj_list)
AdjacencyListToGraph <- function(adj_list){
	nodes <- strsplit(adj_list, " -> ")
	node1 <- sapply(nodes, function(x) x[1])
	node2 <- sapply(nodes, function(x) strsplit(x[2], ","))
	edges <- sapply(seq_along(node1), function(i) cbind(node1[i], node2[[i]])) 
	edges <- do.call(rbind, edges)
	colnames(edges) <- c("node1", "node2")
	nodes <- unique(sort((edges)))
	graph <- list(nodes = nodes, edges = edges)
	return(graph)
}

#' Find all nodes adjacent to a particular node in a graph
#' 
#' \code{ConnectedNodes} finds all the nodes adjacent to \code{node} in the graph with edges specified by code{edges}.
#' 
#' @param node A string giving the name of the node.
#' @param edges A character matrix where each row corresponds to an edge in the graph. It must contain named columns "node1" and "node2" for the outgoing and incoming nodes of the edge (e.g. \code{graph\$edges} for the graph resulting from \code{\link{OverlapGraph}} or \code{\link{DeBruijnGraph}}).
#' @return A character vector listing the nodes connected to \code{node}.
#' @examples
ConnectedNodes <- function(node, edges){
	edge_indices <- which(edges[, "node1"] == node)
	connected_nodes <- edges[edge_indices, "node2"] 
	return(connected_nodes)
}

#' Find an Eulerian cycle in a graph
#' 
#' \code{EulerianCycle} uses the algorithm set forth in the original proof of Euler's Theorem to find an Eulerian cycle in a balanced and strongly connected graph. A graph is balanced if the indegree of every node is equal to its outdegree (see \code{\link{GraphDegree}}) and it is strongly connected if it is possible to reach any node from every other node. An Eulerian cycle is a cycle (i.e. a path that starts and ends at the same node) that traverses each edge of a graph exactly once.
#' 
#' @param graph A (directed) graph in list form as generated by e.g.\code{\link{OverlapGraph}} or \code{\link{AdjacencyListToGraph}}. The graph consists of two named elements, a character vector \code{nodes} containing the nodes and a character matrix \code{edges} containing the edges. The matrix \code{graph\$edges} must contain named columns "node1" and "node2" for the outgoing and incoming nodes of the edges.
#' @return A permutation of \code{graph\$edges} (i.e. a character matrix) that lists the edges in the order in which they are traversed in an Eulerian cycle.
#' @examples
#' adj_list <- c("0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6") 
#' graph <- AdjacencyListToGraph(adj_list)
#' EulerianCycle(graph)
EulerianCycle <- function(graph){
	edges <- graph$edges
	cycle <- NULL
	i <- 1
	repeat {
		# form a cycle until out of options
		while (!is.na(i)) {
			next_node <- edges[i, "node2"]
			cycle <- rbind(cycle, edges[i, ])
			edges <- edges[-i, , drop = FALSE]
			i <- which(edges[, "node1"] == next_node)[1]
		}
		# check if all edges used
		if (nrow(edges) == 0) { 
			break
		}
		# find a node in cycle that has unused edge
		i <- which(edges[, "node1"] %in% unique(as.character(cycle)))[1]
		# check if no remaining node connects to the cycle
		if (is.na(i)) {
			stop("graph is not strongly connected")
		}
		# shift cycle to end at node with unused edge
		next_node <- edges[i, "node1"]
		if (next_node != cycle[1, "node1"]) {
			shift <- which(cycle[, "node1"] == next_node)[1]
			n <- nrow(cycle)
			cycle <- cycle[c(shift:n, 1:(shift - 1)), ]
		}
	}
	return(cycle)
}

#' Print a path in a graph for input into Rosalind
#' 
#' \code{PrintPath} takes a sequence of edges specified by a character matrix which corresponds to a path in a graph and returns a single string specifying the path in the format required for input into Rosalind. For each row in \code{path} the value of "node2" should be equal to the value of "node1" for the subsequent row.
#' 
#' @param path A character matrix which lists the edges of a path. The matrix must contain named columns "node1" and "node2" for the outgoing and incoming nodes of the edges.
#' @return A string which gives the path by listing the nodes of the path separated by " -> ".
#' @examples
#' node1 <- 0:9
#' node2 <- c("3", "0", "1,6", "2", "2", "4", "5,8", "9", "7", "6") 
#' adj_list <- paste(node1, "->", node2)
#' graph <- AdjacencyListToGraph(adj_list)
#' cycle <- EulerianCycle(graph)
#' PrintPath(cycle)
#' 
#' node1 <- c(0:3, 6:9)
#' node2 <- c("2", "3", "1", "0,4", "3,7", "8", "9", "6")
#' adj_list <- paste(node1, "->", node2)
#' graph <- AdjacencyListToGraph(adj_list)
#' path <- EulerianPath(graph)
#' PrintPath(path)
PrintPath <- function(path){
	n <- nrow(path)
	display <- paste0(c(path[, "node1"], path[n, "node2"]), collapse = "->")
	return(display)
}
