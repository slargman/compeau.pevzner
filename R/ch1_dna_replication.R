#' Count the number of occurrences of a given k-mer pattern in a string
#' 
#' \code{PatternCount} returns the number of occurrences of \code{pattern} in the string \code{text}.
#' 
#' @param text character string to be scanned for occurrences of the string \code{pattern}
#' @param pattern character string to be matched in the given character string
#' @examples
#' PatternCount("CGATATATCCATAG", "ATA")
PatternCount <- function(text, pattern){
	count <- 0
	for(i in 1:(nchar(text) - nchar(pattern) + 1)){
		if (identical(substr(text, i, i + nchar(pattern) - 1), pattern)){
			count <- count + 1
		}
	}
	return(count)
}

#' Find the most frequent \emph{k}-mer in a character string
#' 
#' \code{FrequentWords} returns the string of length \emph{k} (\emph{k}-mer) that is the most frequent \emph{k}-mer in the string \code{text}
#' 
#' @param text character string
#' @param k integer which gives the length of \emph{k}-mer
#' @return character vector of all the most frequent \emph{k}-mers in \code{text}
#' @examples
#' FrequentWords("ACAACTATGCATACTATCGGGAACTATCCT", 5)
#' FrequentWords("CGATATATCCATAG", 3)
FrequentWords <- function(text, k){
	count <- vector()
	for (i in 1:(nchar(text) - k + 1)){
		pattern <- substr(text, i, i + k - 1)
		count[i] <- PatternCount(text, pattern)
	}
	max_count <- max(count)
	which_max_count <- which(count == max_count)
	frequent_patterns <- substring(text, which_max_count, which_max_count + k - 1)
	# remove duplicates
	frequent_patterns <- sort(unique(frequent_patterns))
	return(frequent_patterns)
}

#' Transforms nucleotide symbols into lexicographic index
#' 
#' \code{SymbolToNumber} transforms the symbols A, C, G, and T into the respective integers 0, 1, 2, and 3. This corresponds to the indexing starting from 0 used in the book. Inverse of \code{\link{NumberToSymbol}}.
#' 
#' @param symbol character vector with elements equal to either A, C, G, or T
#' @examples
#' SymbolToNumber("A")
#' SymbolToNumber(c("A", "C", "G", "T"))
SymbolToNumber <- function(symbol){
	nucleotides <- c("A", "C", "G", "T")
	if(all(symbol %in% nucleotides)){
		# book uses indexing starting at 0
		number <- (match(symbol, nucleotides) - 1)
		return(number)
	} else {
		return(NA)
	}

}

#' Transforms nucleotide sequence into integer
#' 
#' \code{PatternToNumber} transforms a character string consisting of the symbols A, C, G, and T (eg a DNA sequence) into an integer representing its lexicographic order. Note that the indexing starts at 0 to correspond with the indexing used in the book. Inverse function of \code{\link{NumberToPattern}}.
#' 
#' @param pattern character string consisting only of the characters A, C, G, or T
#' @examples
#' PatternToNumber("AGT")
PatternToNumber <- function(pattern){
	nucleotides <- c("A", "C", "G", "T")
	#check if pattern is empty to establish base case
	if(identical(pattern, "")){
		return(0)
	}
	symbol <- substr(pattern, nchar(pattern), nchar(pattern))
	prefix <- substr(pattern, 1, nchar(pattern) - 1)
	return(4 * PatternToNumber(prefix) + SymbolToNumber(symbol))
}

#' Transforms lexicographic index into corresponding nucleotide symbol
#' 
#' \code{NumberToSymbol} transforms the integers 0, 1, 2, and 3 into the symbols A, C, G, and T with the corresponding lexicographic index. This corresponds to the indexing starting from 0 used in the book. Inverse of \code{\link{SymbolToNumber}}.
#' 
#' @param index integer giving the lexicographic index of the symbol to be returned
#' @examples
#' NumberToSymbol(0)
NumberToSymbol <- function(index){
	nucleotides <- c("A", "C", "G", "T")
	# book uses indexing starting at 0
	return(nucleotides[index + 1])
}

#' Converts integer index to corresponding nucleotide sequence
#' 
#' \code{NumberToPattern} takes an integer \code{index} and converts it to the nucleotide of length \code{k} that has that index among \emph{k}-mers arranged by lexicographic order. Note that indexing starts at 0 to correspond with the indexing convention used in the book. Inverse function of \code{\link{PatternToNumber}}.
#' 
#' @param index integer which give the lexicographic order of the \emph{k}-mer pattern to be found
#' @param k integer which gives length of nucleotide sequence
#' @examples
#' NumberToPattern(11, 3)
NumberToPattern <- function(index, k){
	if(identical(k,1)){
		return(NumberToSymbol(index))
	}
	prefix_index <- index %/% 4
	r <- index %% 4
	symbol <- NumberToSymbol(r)
	prefix_pattern <- NumberToPattern(prefix_index, k - 1)
	return(paste0(prefix_pattern, symbol))
}

#' Generate frequency array giving occurence of all k-mers in text
#' 
#' \code{ComputingFrequencies} takes a character string \code{text} and returns a "frequency array" (numeric vector) which gives the number of occurences of each k-mer in the text. The index of the frequency array corresponds to 1 more than the lexicographic index of the corresponding k-mer (lexicographic index shifted by 1 to correspond with indexing starting at 0 as in the book).
#' 
#' @param text character string to be scanned for frequent k-mers
#' @param k integer which gives length to be used for sequences
#' @examples
#' ComputingFrequencies("AAGCAAAGGTGGG", 2)
ComputingFrequencies <- function(text, k){
	frequency_array <- numeric(4^k)
	for(i in 1:(nchar(text) - k + 1)){
		pattern <- substr(text, i, i + k - 1)
		# shift to account for indexing from 0
		j <- PatternToNumber(pattern) + 1
		frequency_array[j] <- frequency_array[j] + 1
	}
	return(frequency_array)
}

#' Find the most frequent \emph{k}-mer in a character string
#' 
#' \code{FasterFrequentWords} returns the string of length \emph{k} (\emph{k}-mer) that is the most frequent \emph{k}-mer in the string \code{text}. Differs from \code{\link{FrequentWords}} by using a frequency array generated from \code{\link{ComputingFrequencies}}.
#' 
#' @param text character string
#' @param k integer which gives the length of \emph{k}-mer
#' @return character vector of all the most frequent \emph{k}-mers in \code{text}
#' @examples
#' FasterFrequentWords("ACAACTATGCATACTATCGGGAACTATCCT", 5)
#' FasterFrequentWords("CGATATATCCATAG", 3)
#' FasterFrequentWords("AAGCAAAGGTGGG", 2)
FasterFrequentWords <- function(text, k){
	#frequent_patterns <- character(0)
	frequency_array <- ComputingFrequencies(text, k)
	max_count <- max(frequency_array)
	# adjust index by 1
	frequent_patterns <- NumberToPattern(which(frequency_array == max_count) - 1, k)
	return(frequent_patterns)
}

#' Find the most frequent \emph{k}-mer in a character string
#' 
#' \code{FindingFrequentWordsBySorting} returns the string of length \emph{k} (\emph{k}-mer) that is the most frequent \emph{k}-mer in the string \code{text}. Differs from \code{\link{FrequentWords}} and \code{\link{ComputingFrequencies}} by using an index that lists the k-mers that occur in \code{text} in the order they appear using their lexicographic indices. This index is then sorted to find the most frequent \emph{k}-mers.
#' 
#' @param text character string
#' @param k integer which gives the length of \emph{k}-mer
#' @return character vector of all the most frequent \emph{k}-mers in \code{text}
#' @examples
#' FindingFrequentWordsBySorting("ACAACTATGCATACTATCGGGAACTATCCT", 5)
#' FindingFrequentWordsBySorting("CGATATATCCATAG", 3)
#' FindingFrequentWordsBySorting("AAGCAAAGGTGGG", 2)
FindingFrequentWordsBySorting <- function(text, k){
	index <- numeric(0)
	count <- numeric(0)
	for(i in 1:(nchar(text) - k + 1)){
		pattern <- substr(text, i, i + k - 1)
		# list of k-mers that appear in text in order given by lexicographic index of k-mer
		index[i] <- PatternToNumber(pattern)
		count[i] <- 1
	}
	sorted_index <- sort(index)
	for(i in 2:(nchar(text) - k + 1)){
		if(sorted_index[i] == sorted_index[i - 1]){
			count[i] <- count[i - 1] + 1
		}
	}
	max_count <- max(count)
	frequent_patterns <- NumberToPattern(sorted_index[which(count == max_count)], k)
	return(frequent_patterns)
}

#' Give the complement nucleotide for a given nucleotide
#' 
#' \code{NucleotideComplement) returns the nucleotide that is the complement to a given nucleotide, i.e. returns T, G, C, and A for A, C, G, and T respectively. Can operate on a vector, returning a vector of the complements.
#' 
#' @param nucleotide character vector with each element equal to A, C, G, or T
#' @examples
#' NucleotideComplement("A")
#' NucleotideComplement(c("A", "C", "G", "T"))
NucleotideComplement <- function(nucleotide){
	nucleotide_complements <- c("T", "G", "C", "A")
	index <- SymbolToNumber(nucleotide)
	# index starts at 0 for SymbolToNumber
	complement <- nucleotide_complements[index + 1]
	return(complement)
}

#' Give the complement nucleotide sequence for given nucleotide pattern
#' 
#' \code{NucleotideComplement) returns the nucleotide squence that is the reverse complement to a given nucleotide pattern. If the nucleotide pattern is written from 5' to 3', then \code{NucleotideComplement} returns the complementary sequence also written from 5' to 3'. Relies on \code{\link{NucleotideComplement}}.
#' 
#' @param pattern character string consisting of only the characters A, C, G, and T
#' @examples
#' ReverseComplement("AGTCGCATAGT")
ReverseComplement <- function(pattern){
	nucleotides <- strsplit(pattern, "")[[1]]
	nucleotides_reverse_complement <- rev(NucleotideComplement(nucleotides))
	# collapse argument used to collapse vector argument
	pattern_complement <- paste0(nucleotides_reverse_complement, collapse = "")
	return(pattern_complement)
}

#' Find all occurrences of a pattern in a string
#' 
#' \code{PatternMatch} returns the indices (starting at 0) of all starting positions in the string \code{genome} where the string \code{pattern} appears as a substring. Indexing from 0 is used to match the book. Note that code\{PatternMatch} returns overlapping matches unlike \code{\link{gregexpr}} which only returns disjoint matches.
#' 
#' @param pattern string to be searched for
#' @param genome string where occurences of \code{pattern} will be searched for
#' @examples
#' PatternMatch("AA", "CAAAT")
#' PatternMatch("TCA", "ATGATCAAG")
PatternMatch <- function(pattern, genome){
	# gregexpr returns disjoint matches only so perl regular expression with look-ahead assertion (?=...) is necessary
	match_location <- as.integer(gregexpr(paste0("(?=", pattern, ")"), genome, ignore.case = TRUE, perl = TRUE)[[1]])
	# to match indexing from 0 in book
	return(match_location - 1)
}

#' Find patterns forming clumps in a string
#' 
#' \code{WorseClumpFinding} finds all distinct \emph{k}-mers (strings of length \emph{k}) that form (\emph{L},\emph{t})-clumps in \code{genome}. A \emph{k}-mer \code{pattern} forms an (\emph{L},\emph{t})-clump inside a (longer) string \code{genome} if there is an interval of \code{genome} of length \emph{L} in which this \emph{k}-mer appears at least \emph{t} times (assumes that the \emph{k}-mer completely fits in the interval). Differs from \code{\link{ClumpFinding}} in that \code{ClumpFinding} is faster since it only calculates the frequency array using \code{\link{ComputingFrequencies}} once then updates it for each window of \code{genome} of length \emph{L}, whereas \code{WorseClumpFinding} uses \code{\link{ComputingFrequencies}} to calculate the frequency array from scratch for each window of length L.
#' 
#' @param genome character string consisting only of the characters A, C, G, and T
#' @param k integer giving length of k-mers to be checked for clumps
#' @param t integer giving number of occurences of a k-mer in window of length L to qualify as a clump
#' @param L integer giving size of window to be checked for occurences of k-mer to determine presence of clump
#' @examples
#' WorseClumpFinding("GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC", 4, 3, 25)
WorseClumpFinding <- function(genome, k, t, L){
	clump <- logical(4^k)
	# slide window of length L down genome and check for t clumps
	for(i in 1:(nchar(genome) - L + 1)){
		text <- substr(genome, i, i + L - 1)
		frequency_array <- ComputingFrequencies(text, k)
		clump[which(frequency_array >= t)] <- TRUE
	}
	frequent_patterns <- NumberToPattern(which(clump) - 1, k)
	return(frequent_patterns)
}

#' Find patterns forming clumps in a string
#' 
#' \code{ClumpFinding} finds all distinct \emph{k}-mers (strings of length \emph{k}) that form (\emph{L},\emph{t})-clumps in \code{genome}. A \emph{k}-mer \code{pattern} forms an (\emph{L},\emph{t})-clump inside a (longer) string \code{genome} if there is an interval of \code{genome} of length \emph{L} in which this \emph{k}-mer appears at least \emph{t} times (assumes that the \emph{k}-mer completely fits in the interval).
#' 
#' @param genome character string consisting only of the characters A, C, G, and T
#' @param k integer giving length of k-mers to be checked for clumps
#' @param t integer giving number of occurences of a k-mer in window of length L to qualify as a clump
#' @param L integer giving size of window to be checked for occurences of k-mer to determine presence of clump
#' @examples
#' ClumpFinding("GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTTGTTTTAC", 4, 3, 25)
ClumpFinding <- function(genome, k, t, L){
	text <- substr(genome, 1, L)
	frequency_array <- ComputingFrequencies(text, k)
	clump <- (frequency_array >= t)

	# slide window of length L down genome and update frequency_array
	for(i in 2:(nchar(genome) - L + 1)){
		first_pattern <- substr(genome, i - 1, i + k - 2)
		index <- PatternToNumber(first_pattern) + 1
		frequency_array[index] <- frequency_array[index] - 1

		last_pattern <- substr(genome, i + L - k, i + L - 1)
		index <- PatternToNumber(last_pattern) + 1
		frequency_array[index] <- frequency_array[index] + 1
		if(frequency_array[index] >= t){
			clump[index] <- TRUE
		}
	}

	frequent_patterns <- NumberToPattern(which(clump) - 1, k)
	return(frequent_patterns)
}

#' Find location of minimum skew
#' 
#' \code{FindMinimumSkew} finds the location \emph{i} in a character string which has the lowest value of \eqn{\text{Skew}_i(\text{Genome})}. Skew_i(Genome) is defined as the difference between the total number of occurrences of G and the total number of occurrences of C in the first \emph{i} nucleotides of \code{Genome}. Skew_0(Genome) is defined as 0. Note that indexing starts at 1, i.e. FindMiniumSkew("CG") returns 1.
#' 
#' @param genome character string consisting only of the characters A, C, G, or T
#' @examples
#' FindMinimumSkew("CATGGGCATCGGCCATACGCC")
#' FindMinimumSkew("CG")
FindMinimumSkew <- function(genome){
	skew <- numeric(nchar(genome))
	genome_vector <- strsplit(genome, "")[[1]]

	# calculate skew
	for (i in 2:length(genome_vector)) {
		if (genome_vector[i] == "G") {
			skew[i] <- skew[i - 1] + 1
		} else if (genome_vector[i] == "C") {
			skew[i] <- skew[i - 1] - 1
		} else {
			skew[i] <- skew[i - 1]
		}
	}

	min_skew <- min(skew)
	oriC <- which(skew == min_skew)
	return(oriC)
}

#' Compute the Hamming distance between two strings.
#' 
#' \code{HammingDistance} returns the Hamming distance between two string \code{p} and \code{q} of equal length.
#' 
#' We say that the position \emph{i} in \emph{k}-mers \eqn{p_1 \cdots p_k} is a mismatch if \eqn{p_i \neq q_i}. The Hamming distance between \emph{p} and \emph{k} is defined number of mismatches between these two strings.
#' 
#' @param p character string of equal length to q
#' @param q character string of equal length to p
#' @examples
#' HammingDistance("ATGATCAAG", "ATGATCAAC")
#' HammingDistance("ATGATCAAG", "ATGCTCAAC")
#' HammingDistance("ATGATCAAG", "ATGCTCGAC")
HammingDistance <- function(p, q){
	p_vec <- strsplit(p, "")[[1]]
	q_vec <- strsplit(q, "")[[1]]
	mismatches <- sum(p_vec != q_vec)
	return(mismatches)
}

#' Find all approximate occurrences of a pattern in a string
#' 
#' \code{ApproximatePatternMatch} returns the indices (starting at 0) of all starting positions in the string \code{text} where the string \code{pattern} appears as a substring with up to \emph{d} mismatches (i.e. differs \code{pattern} has a maximum Hamming distance of \emph{d} from the substring of \code{text} at that position). Indexing from 0 is used to match the book. Note that code\{ApproximatePatternMatch} returns overlapping matches like \link\code{PatternMatch} and unlike \code{\link{gregexpr}} which only returns disjoint matches.
#' 
#' @param pattern string to be approximately searched for
#' @param text string where occurences of \code{pattern} will be searched for
#' @param d integer giving allowed number of mismatches (i.e. maximum Hamming distance for matches.
#' @examples
#' ApproximatePatternMatch("AAAAA", "AACAAGCATAAACATTAAAGAG", 1)
#' ApproximatePatternMatch("TCA", "ATGATCAAG", 1)
ApproximatePatternMatch <- function(pattern, text, d){
	indices <- 1:(nchar(text) - nchar(pattern) + 1)
	substrings <- substring(text, indices, indices + nchar(pattern) - 1)
	hamming_distance <- numeric(length(substrings))

	# compute Hamming distances
	for (i in 1:length(hamming_distance)){
		hamming_distance[i] <- HammingDistance(pattern, substrings[i])
	}

	matches <- which(hamming_distance <= d)
	# index from 0
	return(matches - 1)
}

#' Count the number of approximate occurrences of a given k-mer pattern in a string
#' 
#' \code{ApproximatePatternCount} returns the number of approximate occurrences of \code{pattern} in the string \code{text}. By approximate occurrences we mean substrings of \code{text} such that the Hamming distance between the substring and \code{pattern} is at most \emph{d} (i.e. there are at most \emph{d} mismatches between the two strings).
#' 
#' @param text character string to be scanned for approximate occurrences of the string \code{pattern}
#' @param pattern character string to be approximately matched in the given character string
#' @param d integer giving the maximum number of allowed mismatches for an approximate match
#' @examples
#' ApproximatePatternCount("AACAAGCATAAACATTAAAGAG", "AAAAA", 1)
ApproximatePatternCount <- function(text, pattern, d){
	count <- 0
	for (i in 1:(nchar(text) - nchar(pattern) + 1)) {
		subtext <- substring(text, i, i + nchar(pattern) - 1)
		if (HammingDistance(pattern, subtext) <= d) {
			count <- count + 1
		}
	}
	return(count)
}

#' Generate the 1-neighborhood of a string
#' 
#' \code{ImmediateNeighbors} returns the set of all strings that are in the 1-neighborhood of the string \code{pattern}, i.e. that differ from \code{pattern} in at most one position. The \emph{d}-neighborhood of a \emph{k}-mer \code{pattern} is the set of all \emph{k}-mers whose Hamming distance from \code{pattern} does not exceed \emph{d}.
#' 
#' This function is used in the definition of \code{\link{IterativeNeighbors}}.
#' 
#' @param pattern A character string.
#' @examples
#' ImmediateNeighbors("ACG")
ImmediateNeighbors <- function(pattern){
	nucleotides <- c("A", "C", "G", "T")
	neighborhood <- pattern
	for (i in 1:nchar(pattern)) {
		symbol <- substring(pattern, i, i)
		# iterate for each nucleotide different from symbol
		for (nuc in nucleotides[nucleotides != symbol]) {
			neighbor <- pattern
			substring(neighbor, i, i) <- nuc
			neighborhood <- c(neighborhood, neighbor)
		}
	}
	return(neighborhood)
}

#' Generate the d-neighborhood of a string
#' 
#' \code{Neighbors} returns the set of all strings that are in the \emph{d}-neighborhood of the string \code{pattern}. The \emph{d}-neighborhood of a \emph{k}-mer \code{pattern} is the set of all \emph{k}-mers whose Hamming distance from \code{pattern} does not exceed \emph{d}.
#' 
#' @param pattern A string.
#' @param d A positive integer giving the maximum Hamming distance for the neighborhood.
#' @examples
#' Neighbors("A", 1)
#' Neighbors("ACG", 2)
Neighbors <- function(pattern, d){
	nucleotides <- c("A", "C", "G", "T")
	# for exact matches
	if (identical(d, 0)){
		return(pattern)
	}
	# nchar returns integer rather than numeric class
	if (identical(nchar(pattern), as.integer(1))){
		return(nucleotides)
	}
	neighborhood <- character(0)
	first_symbol <- substring(pattern, 1, 1)
	suffix <- substring(pattern, 2, nchar(pattern))
	suffix_neighbors <- Neighbors(suffix, d)
	# modify neighborhood of pattern's suffix
	for (text in suffix_neighbors){
		if (HammingDistance(suffix, text) < d) {
			neighborhood <- c(neighborhood, paste0(nucleotides, text))
		} else {
			neighborhood <- c(neighborhood, paste0(first_symbol, text))
		}
	}
	return(neighborhood)
}

#' Generate the d-neighborhood of a string using an iterative algorithm
#' 
#' \code{IterativeNeighbors} returns the set of all strings that are in the \emph{d}-neighborhood of the string \code{pattern}. The \emph{d}-neighborhood of a \emph{k}-mer \code{pattern} is the set of all \emph{k}-mers whose Hamming distance from \code{pattern} does not exceed \emph{d}.
#' 
#' \code{IterativeNeighbors} differs from \code{\link{Neighbors}} in that \code{IterativeNeighbors} uses an iterative algorithm whereas \code{\link{Neighbors}} uses a recursive algorithm.
#' 
#' @param pattern A string.
#' @param d A positive integer giving the maximum Hamming distance for the neighborhood.
#' @examples
#' IterativeNeighbors("A", 1)
#' IterativeNeighbors("ACG", 2)
IterativeNeighbors <- function(pattern, d){
	neighborhood <- pattern
	for (j in 1:d) {
		for (neighbor in neighborhood) {
			neighborhood <- c(neighborhood, ImmediateNeighbors(neighbor))
			neighborhood <- unique(neighborhood)
		}
	}
	return(neighborhood)
}
