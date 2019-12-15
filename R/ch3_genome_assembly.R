#' Generate the k-mer composition of a string.
#' 
#' Given a string \code{text}, its \emph{k}-mer composition is the collection of all \emph{k}-mer substrings of \code{text} (including repeated \emph{k}-mers. \code{StringComposition generates the \emph{k}-mer composition of a given string, listed in lexicographic order.
#' 
#' @param text A string.
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
