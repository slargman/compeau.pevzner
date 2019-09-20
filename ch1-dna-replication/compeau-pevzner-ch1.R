#' Count the number of occurrences of a given k-mer pattern in a string
#' 
#' \code{PatternCount} returns the number of occurrences of \code{pattern} in the string \code{text}
#' 
#' @param text character string to be scanned for occurrences of the string \code{pattern}
#' @param pattern character string to be matched in the given character string
#' @examples
#' PatternCount('CGATATATCCATAG', 'ATA')
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
#' FrequentWords('ACAACTATGCATACTATCGGGAACTATCCT', 5)
#' FrequentWords('CGATATATCCATAG', 3)
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
