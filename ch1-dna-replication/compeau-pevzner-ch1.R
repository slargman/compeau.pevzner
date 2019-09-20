#' Count the number of occurrences of a given k-mer pattern in a string
#' 
#' \code{PatternCount} returns the number of occurrences of \code{pattern} in the string \code{text}
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
#' @param symbol character equal to either A, C, G, or T
#' @examples
#' SymbolToNumber("A")
SymbolToNumber <- function(symbol){
	nucleotides <- c("A", "C", "G", "T")
	if(symbol %in% nucleotides){
		# book uses indexing starting at 0
		number <- (which(symbol == nucleotides) - 1)
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
	#check if pattern is empty
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
