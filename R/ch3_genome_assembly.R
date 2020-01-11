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
#' \code{PathGraph} generates the \code{k}mer path graph for a string \code{text} and specified \code{k}. Given a string \code{text}, the path graph of \code{k}mers for it is the path consisting of \code{nchar(text)} - \code{k} + 1 edges, where the \emph{i}-th edge of this path is labeled by the \emph{i}-th \code{k}-mer in \code{text} and the \emph{i}-th node of the path is labeled by the \emph{i}-th (\code{k} - 1)-mer in \code{text}.
#' 
#' @inheritParams StringComposition
#' @return A list representing the path graph of the collection of \emph{k}-mers \code{pattern}. This list has two named elements. The element \code{graph\$nodes} contains a character vector listing the nodes of the graph. The element \code{graph\$edges} contains a three column character matrix where each row corresponds to an edge in the graph and the columns list respectively the outgoing node, the edge label, and the incoming node of an edge.
#' @examples
#' text <- "TAATGCCATGGGATGTT"
#' k <- 3
#' PathGraph(text, k)
PathGraph <- function(text, k){
	n <- nchar(text)
	# edges labelled by kmers
	edge_labels <- substring(text, 1:(n - k + 1), k:n)
	# nodes labelled by (k - 1)mers
	nodes <- substring(text, 1:(n - k + 2), (k - 1):n)
	m <- length(nodes)
	edges <- matrix(c(nodes[1:(m - 1)], edge_labels, nodes[2:m]), ncol = 3)
	colnames(edges) <- c("node1", "label", "node2")
	lexicographic <- order(edges[, "node1"], edges[, "node2"])
	edges <- edges[lexicographic, ]
	graph <- list(nodes = nodes, edges = edges)
	return(graph)
}
