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
