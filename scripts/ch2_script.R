
dataset <- readLines("../data/rosalind_ba2a.txt")
k <- as.integer(substring(dataset[1], 1, 1)[[1]])
d <- as.integer(substring(dataset[1], 3, 3)[[1]])
dna <- dataset[-1]
answer <- MotifEnumeration(dna, k, d)
output <- paste(answer, collapse = " ")
writeClipboard(output)
