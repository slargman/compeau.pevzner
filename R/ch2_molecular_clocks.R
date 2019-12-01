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
