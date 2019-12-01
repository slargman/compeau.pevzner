#' Find all (\emph{k}, \emph{d})-motifs in a collection of strings
#' 
#' \code{MotifEnumeration} finds all (\emph{k}, \emph{d})-motifs in a collection of strings \code{dna}. Given a collection of strings \code{dna} and an integer \emph{d}, a \emph{k}-mer is a (\emph{k}, \emph{d})-motif if it appears in every string from \code{dna} with at most \emph{d} mismatches.
#' 
#' @param dna A character vector containing strings to be searched for motifs.
#' @param k An integer which gives the length of motifs.
#' @param d An integer which gives the number of mismatches allowed in a motif.
#' @return A character vector containing all (\emph{k},\emph{d})-motifs in \code{dna}.
#' @examples
#' MotifEnumeration(c("ATCGA", "CCAGC", "GCATG"), 3, 1)
#' MotifEnumeration(c("ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"), 3, 1)
MotifEnumeration <- function(dna, k, d){
	patterns <- character(0)
	for (string in dna) {
		# for each k-mer pattern in dna
		for (i in 1:(nchar(string) - k + 1)) {
			pattern <- substring(string, i, i + k - 1)
			neighborhood <- Neighbors(pattern, d) 
			# for each d-neighbor of pattern
			for (neighbor in neighborhood) {
				neighbor_count <- integer(0)
				for (j in seq_along(dna)) {
					neighbor_count[j] <- ApproximatePatternCount(dna[j], neighbor, d)
				}
				# if pattern appears in each string of dna with at most d mismatches
				if (all(neighbor_count > 0)) {
					patterns <- c(patterns, neighbor)
				}
			}
		}
	}
	patterns <- unique(patterns)
	return(patterns)
}

#' Compute the distance between a string and a collection of potentially longer strings
#' 
#' \code{StringDistance} computes the distance between two strings, which may be of different lengths, or between a string and a collection of strings of equal or greater length. If \code{p} and \code{q} are strings of equal length, the distance \emph{d}(\emph{p}, \emph{q}) is given by the Hamming distance (implemented in \code{\link{HammingDistance}}). If \code{p} is a \emph{k}-mer and \code{q} is a longer string than \code{p}, then \emph{d}(\emph{p}, \emph{q})) is defined as the minimum Hamming distance between \emph{p} and any \emph{k}-mer in \emph{q}. Finally, given a \emph{k}-mer \emph{p} and a set of strings \emph{q} (i.e. a character vector), we define \emph{d}(\emph{p}, \emph{q}) as the sum of distances between \emph{p} and all strings in \emph{q}.
#' 
#' Only one of \code{p} and \code{q} can have length greater than 1. If \code{p} has length greater than 1, each element must have at least as many characters as \code{q}, and vice versa.
#' 
#' @param p A character vector.
#' @param q A character vector.
#' @return An integer giving the distance between \code{p} and \code{q}.
#' @examples
#' StringDistance("ATGATCAAG", "ATGATCAAC")
#' StringDistance("ATGATCAAG", "ATGCTCAAC")
#' StringDistance("ATGATCAAG", "ATGCTCGAC")
#' StringDistance("GATTCTCA", "GCAAAGACGCTGACCAA")
#' StringDistance("GCAAAGACGCTGACCAA", "GATTCTCA")
#' StringDistance("GATTCTCA", c("GCAAAGACGCTGACCAA", "TCGAACGACTCAGTCGATCA"))
#' StringDistance(c("GCAAAGACGCTGACCAA", "TCGAACGACTCAGTCGATCA"), "GATTCTCA")
StringDistance <- function(p, q){
	if (length(p) > 1 && length(q) > 1) {
		stop("only one argument can have length greater than 1")
	} else if (length(p) > 1) {
		pattern <- q
		text <- p
	} else if (length(q) > 1) {
		pattern <- p
		text <- q
	} else if (nchar(p) >= nchar(q)) {
		pattern <- q
		text <- p
	} else if (nchar(q) > nchar(p)) {
		pattern <- p
		text <- q
	}
	if (any(nchar(text) < nchar(pattern))) {
		stop("strings in the collection cannot be shorter than the comparison string")
	}
	distance <- 0
	k <- nchar(pattern)
	for (i in 1:length(text)) {
		i_distance <- Inf
		for (j in 1:(nchar(text[i]) - k + 1)) {
			subtext <- substring(text[i], j, j + k - 1)
			subtext_distance <- HammingDistance(pattern, subtext)
			if (subtext_distance < i_distance) {
				i_distance <- subtext_distance
			}
		}
		distance <- distance + i_distance
	}
	return(distance)
}

#' Find a median string for a collection of strings
#' 
#' \code{MedianString} finds the median string of length \emph{k} for the collection of strings (i.e. character vector) \code{dna}. A median string for a collection of strings \code{dna} is a \emph{k}-mer \code{pattern} that minimizes \emph{d}(\code{pattern}, \code{dna}) (see \code{\link{StringDistance}}) over all \emph{k}-mers \code{pattern}.
#' 
#' @param dna A character vector (i.e. a collection of character strings).
#' @param k A positive integer giving the length for the median string.
#' @return A string which is the median string for \code{dna}.
#' @examples
#' MedianString(c("AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTACGGGACAG"), 3)
MedianString <- function(dna, k){
	distance <- Inf
	for (index in 0:(4^k - 1)) {
		pattern <- NumberToPattern(index, k)
		pattern_distance <- StringDistance(pattern, dna)
		if (pattern_distance < distance) {
			distance <- pattern_distance
			median_string <- pattern
		}
	}
	return(median_string)
}
