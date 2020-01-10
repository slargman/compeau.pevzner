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

#' Find a Profile-most probable \emph{k}-mer in a string
#' 
#' \code{ProfileMostProbableString} finds the \emph{k}-mer substring in a string \code{text} that is the most probable string according to the probability distribution defined in \code{profile}. The matrix \code{profile} is 4 by \emph{k}, with the \emph{i}th column giving the probability distribution for the \emph{i}th character in a \emph{k}-mer and the rows corresponding to probabilities for A, C, G, and T respectively.
#' 
#' If there are multiple \code{profile}-most probable \emph{k}-mers in \code{text}, then the first such \emph{k}-mer occurring in \code{text} is selected.
#' 
#' @param text A string to be searched for the profile-most probable substring.
#' @param k An integer giving the length of the profile-most probable substring.
#' @param profile A 4 by \emph{k} matrix giving the probability of each nucleotide in a position.
#' @return A string giving a \code{profile}-most probable \emph{k}-mer in \code{text}.
#' @examples
#' text <- "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
#' k <- 5
#' A <- c(0.2, 0.2, 0.3, 0.2, 0.3)
#' C <- c(0.4, 0.3, 0.1, 0.5, 0.1)
#' G <- c(0.3, 0.3, 0.5, 0.2, 0.4)
#' T <- c(0.1, 0.2, 0.1, 0.1, 0.2)
#' profile <- rbind(A, C, G, T)
#' ProfileMostProbableString(text, k, profile)
ProfileMostProbableString <- function(text, k, profile){
	probability <- numeric(0)
	subtext_probability <- numeric(k)
	for (i in 1:(nchar(text) - k + 1)) {
		subtext <- substring(text, i, i + k - 1)
		subtext_vector <- strsplit(subtext, split = '')[[1]]
		subtext_numbers <- SymbolToNumber(subtext_vector) + 1
		for (j in 1:k) {
			subtext_probability[j] <- profile[subtext_numbers[j], j]
		}
		probability[i] <- prod(subtext_probability)
	}
	max_probability <- max(probability)
	index <- which(probability == max_probability)[1]
	most_probable_string <- substring(text, index, index + k - 1)
	return(most_probable_string)
}

#' Generate the motif matrix for a collection of strings
#' 
#' \code{MotifMatrix} forms the motif matrix for a collection (i.e. character vector) of DNA strings of equal length. Given a vector containing \emph{t} DNA strings each of length \emph{k}, \code{MotifMatrix} returns a \emph{t} by \emph{k} matrix where each row corresponds to a string and each column to a character position for the strings.
#' 
#' @param motifs A character vector where each element has the same number of characters.
#' @return A character matrix corresponding to the motif matrix.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' MotifMatrix(motifs)
MotifMatrix <- function(motifs){
	motif_list <- strsplit(motifs, split = "")	
	motif_matrix <- do.call(rbind, motif_list)
	return(motif_matrix)
}

#' Transform a motif matrix into a collection of strings
#' 
#' \code{MotifString} takes a motif matrix and collapses each row to generate the string represented by that row. It returns a character vector where each element of the vector corresponds to a row of the motif matrix.
#' 
#' If \code{motif_matrix} is \emph{t} by \emph{k}, \code{MotifString} produces a vector of length \emph{t}, where each element has \emph{k} characters.
#' 
#' \code{MotifString} is the inverse of \code{\link{MotifMatrix}}.
#' @param motif_matrix A character matrix containing only A, C, G, or T.
#' @return A character vector.
#' @examples
#' motifs1 <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motif_matrix <- MotifMatrix(motifs1)
#' motifs2 <- MotifString(motif_matrix)
MotifString <- function(motif_matrix){
	motifs <- apply(motif_matrix, MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
	return(motifs)
}

#' Count occurences of each nucleotide in a motif matrix column
#' 
#' \code{ColumnCount} returns a 4 by 1 matrix (i.e. a column) that lists the number of occurences of each nucleotide in a character vector or one column matrix consisting only of the characters A, C, G, or T. The rows of the output correspond respectively to the number of occurences of A, C, G, and T in \code{column}.
#' 
#' \code{ColumnCount} is used to implement \code{\link{MotifCount}}.
#' 
#' @param column A character vector or one column matrix.
#' @return A 4 by 1 matrix giving the number of nucleotide occurences.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' ColumnCount(motifs[, 1])
ColumnCount <- function(column){
	A_count <- sum(column == "A")
	C_count <- sum(column == "C")
	G_count <- sum(column == "G")
	T_count <- sum(column == "T")
	count <- matrix(c(A_count, C_count, G_count, T_count))
	return(count)
}

#' Count occurences of each nucleotide in a motif matrix
#' 
#' \code{MotifCount} returns a 4 row matrix that lists the number of occurences of each nucleotide in a motif matrix, i.e. a character matrix consisting only of the characters A, C, G, or T. If the strings in \code{motifs} are \emph{k}-mers, i.e. \code{motifs} has \emph{k} columns, then \code{MotifCount} returns a 4 by \emph{k} matrix. The rows of the output denote respectively the number of occurences of A, C, G, and T in the corresponding column of \code{motifs}. The (\emph{i}, \emph{j})-th element of \code{MotifCount(motifs)} stores the number of times that nucleotide \emph{i} appears in column \emph{j} of \code{motifs}.
#' 
#' \code{ColumnCount} is used to implement \code{\link{MotifCount}}.
#' 
#' @param motifs A character matrix containing only A, C, G, or T.
#' @return A 4 row matrix giving the number of nucleotide occurences.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' MotifCount(motifs)
MotifCount <- function(motifs){
	# vector returns NULL for nrow()
	if (is.null(dim(motifs))) {
		stop("MotifCount requires a matrix")
	}
	motif_count <- apply(motifs, MARGIN = 2, FUN = ColumnCount)
	row.names(motif_count) <- c("A", "C", "G", "T")
	return(motif_count)
}

#' Determine the consensus string for a motif matrix
#' 
#' \code{MotifConsensus} forms the consensus string from the most popular nucleotides in each column of the motif matrix \code{motifs} (ties are awarded to the first nucleotide by lexicographic order). The most popular nucleotide for a position is the one with the highest number of counts in the count matrix (generated by \code{\link{MotifCount}) in the corresponding column.
#' @param motifs A character matrix containing only A, C, G, or T.
#' @return A string.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' MotifConsensus(motifs[, 1])
#' MotifConsensus(motifs)
MotifConsensus <- function(motifs){
	motifs <- as.matrix(motifs)
	motif_count <- MotifCount(motifs)
	motif_max <- apply(motif_count, 2, function(x) which(x == max(x))[1] - 1)
	consensus <- paste0(NumberToSymbol(motif_max), collapse = "")
	return(consensus)
}

#' Calculate the score for a motif matrix column
#' 
#' \code{ColumnScore} returns the score for the one column motif matrix \code{motifs}. The score for a motif matrix is the number of unpopular letters in the motif matrix, where an unpopular letter is one that differs from the most frequent nucleotide for that position. If there are two most frequent nucleotides the first one is awarded the tie, as in the determination for \code{\link{MotifConsensus}}.
#' 
#' \code{\link{ColumnScore}} is used to implement \code{MotifScore}.
#' 
#' @param column A one column character matrix containing only A, C, G, or T.
#' @return The score for the motif column matrix \code{column}.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' ColumnScore(motifs[, 1])
ColumnScore <- function(column){
	consensus <- MotifConsensus(column)
	column[which(column != consensus)] <- 1
	column[which(column == consensus)] <- 0
	score <- sum(as.integer(column))
	return(score)
}

#' Calculate the score for a motif matrix
#' 
#' \code{MotifScore} returns the score for the motif matrix \code{motifs}. The score for a motif matrix is the number of unpopular letters in the motif matrix, where an unpopular letter is one that differs from the most frequent nucleotide for that position. If there are two most frequent nucleotides the first one is awarded the tie, as in the determination for \code{\link{MotifConsensus}}.
#' 
#' \code{\link{ColumnScore}} is used to implement \code{MotifScore}.
#' 
#' @param motifs A character matrix containing only A, C, G, or T.
#' @return The score for the motif matrix \code{motifs}.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' MotifScore(motifs)
MotifScore <- function(motifs){
	score <- apply(motifs, MARGIN = 2, FUN = ColumnScore)
	score <- sum(score)
	return(score)
}

#' Generate the profile matrix for a motif matrix
#' 
#' \code{MotifProfile} generates the profile matrix for a given motif matrix. The (\emph{i}, \emph{j})th entry of the profile matrix is the frequency (or probability) of the \emph{i}-th nucleotide in the \emph{j}-th column of the motif matrix. Note that the elements of any column of the profile matrix sum to 1. The profile matrix is generated from the count matrix (produced by \code{\link{MotifCount}}) by dividing all the elements in the count matrix by \emph{t}, the number of rows in \code{motifs}.
#' 
#' See \code{\link{MotifProfilePseudcounts}} for a version that uses pseudocounts to ensure that there are no zeros in the profile matrix.
#'
#' @param motifs A character matrix containing only A, C, G, or T.
#' @return A numeric matrix corresponding to the profile matrix for \code{motifs}.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' MotifProfile(motifs)
MotifProfile <- function(motifs){
	# vector returns NULL for nrow()
	if (is.null(dim(motifs))) {
		stop("MotifProfile requires a matrix")
	}
	counts <- MotifCount(motifs)
	rows <- nrow(motifs)
	motif_profile <- counts / rows
	return(motif_profile)
}

#' Generate the profile matrix for a motif matrix using pseudocounts
#' 
#' \code{MotifProfilePseudocounts} generates the profile matrix for a given motif matrix after incorporating pseudocounts. The (\emph{i}, \emph{j})th entry of the profile matrix is the frequency (or probability) of the \emph{i}-th nucleotide in the \emph{j}-th column of the motif matrix. Note that the elements of any column of the profile matrix sum to 1.
#' 
#' \code{MotifProfilePseudocounts} differs from \code{\link{MotifProfile}} in that it incorporates pseudocounts to ensure that the profile matrix does not contain any zeros (i.e. adding 1 to every element of the counts matrix). The profile matrix is generated from the count matrix (produced by adding one to result from \code{\link{MotifCount}}) by dividing all the elements in the count matrix by \emph{t}, the number of rows in \code{motifs}.
#'
#' @param motifs A character matrix containing only A, C, G, or T.
#' @return A numeric matrix corresponding to the profile matrix for \code{motifs} after incorporating pseudcounts.
#' @examples
#' motifs <- c("ATGAC","CGATG","CTGAT","ACGTC", "CACAC")
#' motifs <- MotifMatrix(motifs)
#' MotifProfilePseudocounts(motifs)
MotifProfilePseudocounts <- function(motifs){
	# vector returns NULL for nrow()
	if (is.null(dim(motifs))) {
		stop("MotifProfilePseudocounts requires a matrix")
	}
	# add 1 for pseudocounts
	counts <- MotifCount(motifs) + 1
	# adjust for added pseudocounts
	rows <- nrow(motifs) + 4
	motif_profile <- counts / rows
	return(motif_profile)
}

#' Generate motifs from a collection of strings and a motif profile
#' 
#' Given a collection of strings \code{dna} and a numeric 4 by \emph{k} matrix \code{motif_profile} (a motif profile), \code{Motifs(motif_profile, dna) returns the collection of \emph{k}-mers formed by the profile-most probable \emph{k}-mers in each sequence from \code{dna}. The profile-most probable string is found using \code{\link{ProfileMostProbableString}}.
#' 
#' \code{Motifs} can accept either a collection of strings or the corresponding motif matrix generated by \code{\link{MotifMatrix}} as input for \code{dna}.
#' 
#' @param motif_profile A 4 by k numeric matrix whose columns sum to 1.
#' @param dna A character vector where each element has at least \emph{k} characters or a character matrix with at least \emph{k} columns.
#' @return A character vector containing the profile-most probable strings in each element of \code{dna}.
#' @examples
#' dna <- c("TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT")
#' dna_matrix <- MotifMatrix(dna)
#' A <- c(.8, 0, 0, .2)
#' C <- c(0, .6, .2, 0)
#' G <- c(.2, .2, .8, 0)
#' T <- c(0, .2, 0, .8)
#' motif_profile <- rbind(A, C, G, T)
#' Motifs(motif_profile, dna)
#' Motifs(motif_profile, MotifMatrix(dna))
Motifs <- function(motif_profile, dna){
	# check for motif matrix
	if (!is.null(dim(dna))) {
		dna <- MotifString(dna)
	}
	k <- dim(motif_profile)[2]
	motifs <- sapply(dna, FUN = function(text) ProfileMostProbableString(text, k, motif_profile), USE.NAMES = FALSE)
	motifs <- MotifMatrix(motifs)
	return(motifs)
}

#' Find motifs in a collection of strings using a greedy algorithm
#' 
#' \code{GreedyMotifSearch} searches a collection of strings \code{dna} for a motif \code{k} characters long using a greedy algorithm (see \code{\link{MotifEnumeration}} for an explanation of motifs). \code{GreedyMotifSearch} returns the motif collection searched by the greedy algorithm that has the lowest \link[MotifScore]{score}. The output is represented by its motif matrix (as in \code{\link{MotifMatrix}}) where each row represents a motif occurrence.
#' 
#' @param dna A character vector where each element has at least \code{k} characters.
#' @param k An integer giving the numbers of characters in the motifs.
#' @param t An integer. Specifies that the first \code{t} elements of \code{dna} should be searched for a motif (defaults to all elements).
#' @return A character matrix. Each row gives the occurrence of the motif in the corresponding element of \code{dna}. The output is a motif matrix like that generated by \code{\link{MotifMatrix}}.
#' @examples
#' k <- 3
#' t <- 5
#' dna <- c("GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG")
#' GreedyMotifSearch(dna, k, t)
GreedyMotifSearch <- function(dna, k, t = length(dna)){
	best_motifs <- MotifMatrix(substring(dna, 1, k))
	best_score <- MotifScore(best_motifs)
	for (l in 1:(nchar(dna[1]) - k + 1)) {
		motifs <- MotifMatrix(substring(dna[1], l, l + k - 1))
		for (i in 2:t) {
			motif_profile <- MotifProfile(motifs[1:(i - 1), , drop = FALSE])
			new_motif <- MotifMatrix(ProfileMostProbableString(dna[i], k, motif_profile))
			motifs <- rbind(motifs, new_motif)
		}
		if (MotifScore(motifs) < best_score) {
			best_motifs <- motifs
			best_score <- MotifScore(best_motifs)
		}
	}
	return(best_motifs)
}

#' Find motifs using a greedy algorithm and pseudocounts
#' 
#' \code{GreedyMotifSearchPseudocounts} searches a collection of strings \code{dna} for a motif \code{k} characters long using a greedy algorithm (see \code{\link{MotifEnumeration}} for an explanation of motifs) and pseudocounts. \code{GreedyMotifSearchPseudocounts} returns the motif collection searched by the greedy algorithm that has the lowest \link[MotifScore]{score}. The output is represented by its motif matrix (as in \code{\link{MotifMatrix}}) where each row represents a motif occurrence.
#' 
#' \code{GreedyMotifSearchPseudocounts} differs from \code{\link{GreedyMotifSearch}} in that it incorporates pseudocounts to ensure that the profile matrix does not contain any zeros (i.e. adding 1 to every element of the counts matrix).
#' 
#' @param dna A character vector where each element has at least \code{k} characters.
#' @param k An integer giving the numbers of characters in the motifs.
#' @param t An integer. Specifies that the first \code{t} elements of \code{dna} should be searched for a motif (defaults to all elements).
#' @return A character matrix. Each row gives the occurrence of the motif in the corresponding element of \code{dna}. The output is a motif matrix like that generated by \code{\link{MotifMatrix}}.
#' @examples
#' k <- 3
#' t <- 5
#' dna <- c("GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG")
#' GreedyMotifSearchPseudocounts(dna, k, t)
GreedyMotifSearchPseudocounts <- function(dna, k, t = length(dna)){
	best_motifs <- MotifMatrix(substring(dna, 1, k))
	best_score <- MotifScore(best_motifs)
	for (l in 1:(nchar(dna[1]) - k + 1)) {
		motifs <- MotifMatrix(substring(dna[1], l, l + k - 1))
		for (i in 2:t) {
			motif_profile <- MotifProfilePseudocounts(motifs[1:(i - 1), , drop = FALSE])
			new_motif <- MotifMatrix(ProfileMostProbableString(dna[i], k, motif_profile))
			motifs <- rbind(motifs, new_motif)
		}
		if (MotifScore(motifs) < best_score) {
			best_motifs <- motifs
			best_score <- MotifScore(best_motifs)
		}
	}
	return(best_motifs)
}

#' Randomly select a substring of a particular length
#' 
#' \code{RandomSubstring} randomly selects a substring in each element of a character vector \code{strings}. The selection in each element is independent and a vector \code{prob} of probability weights for the selectable positions for the start of the substrings can be specified.
#' 
#' @param strings A character vector where each element has at least \code{k} characters.
#' @param k An integer or integer vector giving the numbers of characters in the substrings.
#' s1 <- "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
#' s2 <- "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
#' s3 <- "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
#' s4 <- "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
#' s5 <- "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
#' strings <- c(s1, s2, s3, s4, s5)
#' k <- 8
#' nchar(s1)
#' substrings <- RandomSubstring(strings, k)
#' substrings <- RandomSubstring(strings, k, prob = c(1, rep(0, 24)))
RandomSubstring <- function(strings, k, prob = NULL){
	strings_nchar <- nchar(strings) - k + 1
	first <- sapply(strings_nchar, FUN = function(x) sample(x = x, size = 1, prob = prob))
	motifs <- substring(strings, first = first, last = first + k - 1)
}

#' Find motifs using a randomized algorithm
#' 
#' \code{RandomizedMotifSearch} searches a collection of strings \code{dna} for a motif \code{k} characters long using a randomized algorithm (see \code{\link{MotifEnumeration}} for an explanation of motifs) and pseudocounts (see \code{\link{MotifProfilePseudocounts}}). \code{RandomizedMotifSearch} returns a motif collection as represented by its motif matrix (as in \code{\link{MotifMatrix}}) where each row represents the string correpsonding to a motif occurrence.
#'  
#' @param dna A character vector where each element has at least \code{k} characters.
#' @param k An integer giving the numbers of characters in the motifs.
#' @param t An integer. Specifies that the first \code{t} elements of \code{dna} should be searched for a motif (defaults to all elements).
#' @param iter An integer giving the number of iterations of the algorithm to run. The motifs with the best score (see \code{\link{MotifScore}}) over all iterations is returned.
#' @return A character matrix. Each row gives the occurrence of the motif in the corresponding element of \code{dna}. The output is a motif matrix like that generated by \code{\link{MotifMatrix}}.
#' @examples
#' k <- 8
#' t <- 5
#' s1 <- "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
#' s2 <- "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
#' s3 <- "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
#' s4 <- "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
#' s5 <- "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
#' dna <- c(s1, s2, s3, s4, s5)
#' motifs <- RandomizedMotifSearch(dna, k, t, iter = 1000)
RandomizedMotifSearch <- function(dna, k, t = length(dna), iter = 1){
	top_motifs <- MotifMatrix(RandomSubstring(dna, k))
	for (i in 1:iter) {
		# elements of dna might have different lengths
		motifs <- MotifMatrix(RandomSubstring(dna, k))
		best_motifs <- motifs
		repeat {
			motif_profile <- MotifProfilePseudocounts(motifs)
			motifs <- Motifs(motif_profile, dna)
			# continue loop while scoring is improving
			if (MotifScore(motifs) < MotifScore(best_motifs)) {
				best_motifs <- motifs
			} else {
				break
			}
		}
		if (MotifScore(best_motifs) < MotifScore(top_motifs)) {
			top_motifs <- best_motifs
		}
	}
	return(top_motifs)
}

#' Select a substring according to the distribution in a profile matrix
#' 
#' \code{ProfileRandomString} randomly selects a substring of the string \code{dna} according to the probability distribution implied by the profile matrix \code{motif_profile} (see \code{\link{MotifProfile}}). The conditional probability of each substring of \code{dna} is determined using the profile matrix and these conditional probabilities are normalized to a probability mass function which is used to randomly select a substring of \code{dna}.
#' 
#' \code{ProfileRandomString} can accept \code{dna} as either a string or the equivalent one row motif matrix (\code{\link{MotifMatrix}}).
#' 
#' @param dna A string containing at least \code{k} characters or a one row motif (character) matrix with at least \code{k} columns.
#' @param motif_profile A 4 by k numeric matrix whose columns sum to 1.
#' @param k An integer giving the numbers of characters in the substring.
#' @return A string (not the corresponding motif matrix).
#' @examples
#' s1 <- "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
#' s2 <- "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
#' s3 <- "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
#' s4 <- "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
#' s5 <- "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
#' strings <- c(s1, s2, s3, s4, s5)
#' motif_matrix <- MotifMatrix(strings)
#' motif_profile <- MotifProfile(motif_matrix)
#' k <- 8
#' sub1 <- ProfileRandomString(s1, motif_profile, k)
#' sub2 <- ProfileRandomString(motif_matrix[1, , drop = F], motif_profile, k)
ProfileRandomString <- function(dna, motif_profile, k){
	nucleotides <- 1:4
	names(nucleotides) <- c("A", "C", "G", "T")
	#ensure we have a matrix
	if (is.null(dim(dna))) {
		motif_matrix <- MotifMatrix(dna)
	} else {
		motif_matrix <- dna
	}
	motif_string <- MotifString(motif_matrix)
	n <- ncol(motif_matrix) - k + 1
	substring_probability <- numeric(n)
	for(i in 1:n) {
		nucleotides_num <- nucleotides[motif_matrix[1, i:(i + k - 1)]]
		nucleotide_probability <- motif_profile[matrix(c(nucleotides_num, 1:k), ncol = 2)]
		substring_probability[i] <- prod(nucleotide_probability)
	}
	random_string <- RandomSubstring(motif_string, k, prob = substring_probability)
	return(random_string)
}

#' Find motifs using Gibbs sampling
#' 
#' \code{GibbsSampler} searches a collection of strings \code{dna} for a motif \code{k} characters long using a randomized algorithm (see \code{\link{MotifEnumeration}} for an explanation of motifs) and pseudocounts (see \code{\link{MotifProfilePseudocounts}}). \code{GibbsSampler} returns a motif collection as represented by its motif matrix (as in \code{\link{MotifMatrix}}) where each row represents the string correpsonding to a motif occurrence.
#'  
#' @param dna A character vector where each element has at least \code{k} characters.
#' @param k An integer giving the numbers of characters in the motifs.
#' @param t An integer. Specifies that the first \code{t} elements of \code{dna} should be searched for a motif (defaults to all elements).
#' @param iter An integer giving the number of iterations of the algorithm to run. The motifs with the best score (see \code{\link{MotifScore}}) over all iterations is returned.
#' @return A character matrix. Each row gives the occurrence of the motif in the corresponding element of \code{dna}. The output is a motif matrix like that generated by \code{\link{MotifMatrix}}.
#' @examples
#' k <- 8
#' t <- 5
#' N <- 100
#' s1 <- "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA"
#' s2 <- "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG"
#' s3 <- "TAGTACCGAGACCGAAAGAAGTATACAGGCGT"
#' s4 <- "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC"
#' s5 <- "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
#' dna <- c(s1, s2, s3, s4, s5)
#' motifs <- GibbsSampler(dna, k, t, N, iter = 100)
GibbsSampler <- function(dna, k, t = length(dna), N, iter = 1){
	top_motifs <- MotifMatrix(RandomSubstring(dna, k))
	for (r in 1:iter) {
		motifs <- MotifMatrix(RandomSubstring(dna, k))
		best_motifs <- motifs
		for (j in 1:N) {
			i <- sample(t, size = 1)
			motif_profile <- MotifProfilePseudocounts(motifs[-i, , drop = FALSE])
			motifs[i, ] <- MotifMatrix(ProfileRandomString(dna[i], motif_profile, k))
			if (MotifScore(motifs) < MotifScore(best_motifs)) {
				best_motifs <- motifs
			}
		}
		if (MotifScore(best_motifs) < MotifScore(top_motifs)) {
			top_motifs <- best_motifs
		}
	}
	return(top_motifs)
}
