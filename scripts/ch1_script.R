# 1B
dataset <- readLines("../data/rosalind_ba1b.txt")
text <- dataset[[1]]
k <- as.integer(dataset[[2]])
frequent_words <- FrequentWords(text,k)
output <- paste(frequent_words, collapse = " ")
writeClipboard(output)

# 1C
dataset <- readLines("../data/rosalind_ba1c.txt")
pattern <- dataset[[1]]
complement <- ReverseComplement(pattern)
output <- paste(complement, collapse = " ")
writeClipboard(output)

# 1D
dataset <- readLines("../data/rosalind_ba1d.txt")
pattern <- dataset[[1]]
genome <- dataset[[2]]
occurences <- PatternMatch(pattern, genome)
output <- paste(occurences, collapse = " ")
writeClipboard(output)

writeClipboard(paste0(as.character(PatternMatch(readLines("dataset_3_5.txt")[[1]], readLines("dataset_3_5.txt")[[2]])), collapse = " "))

cholerae <- readLines("vibrio_cholerae.txt")
