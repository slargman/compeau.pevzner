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

# 1E
dataset <- readLines("../data/rosalind_ba1e.txt")
genome <- dataset[[1]]
vars <- as.integer(strsplit(dataset[[2]], split = " ")[[1]])
k <- vars[1]
L <- vars[2]
t <- vars[3]
clumps <- ClumpFinding(genome, k, t, L)
output <- paste(clumps, collapse = " ")
writeClipboard(output)


# 1F
dataset <- readLines("../data/rosalind_ba1f.txt")
genome <- dataset[[1]]
answer <- FindMinimumSkew(genome)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1G
dataset <- readLines("../data/rosalind_ba1g.txt")
s1 <- dataset[[1]]
s2 <- dataset[[2]]
answer <- HammingDistance(s1, s2)
output <- paste(answer, collapse = " ")
writeClipboard(output)


# 1H
dataset <- readLines("../data/rosalind_ba1h.txt")
#dataset <- list("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC", 3)
pattern <- dataset[[1]]
text <- dataset[[2]]
d <- as.integer(dataset[[3]])
answer <- ApproximatePatternMatch(pattern, text, d)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1I
dataset <- readLines("../data/rosalind_ba1i.txt")
text <- dataset[[1]]
vars <- as.integer(strsplit(dataset[[2]], split = " ")[[1]])
k <- vars[1]
d <- vars[2]
answer <- FrequentWordsWithMismatches(text, k, d)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1J
dataset <- readLines("../data/rosalind_ba1j.txt")
#dataset <- list("ACGTTGCATGTCGCATGATGCATGAGAGCT", "4 1")
text <- dataset[[1]]
vars <- as.integer(strsplit(dataset[[2]], split = " ")[[1]])
k <- vars[1]
d <- vars[2]
answer <- FrequentWordsWithMismatchesAndComplements(text, k, d)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1K
dataset <- readLines("../data/rosalind_ba1k.txt")
text <- dataset[[1]]
k <- as.integer(dataset[[2]])
answer <- ComputingFrequencies(text, k)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1L
dataset <- readLines("../data/rosalind_ba1l.txt")
#dataset <- "AGT"
pattern <- dataset[[1]]
answer <- PatternToNumber(pattern)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1M
dataset <- readLines("../data/rosalind_ba1m.txt")
index <- as.integer(dataset[[1]])
k <- as.integer(dataset[[2]])
answer <- NumberToPattern(index, k)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 1N
dataset <- readLines("../data/rosalind_ba1n.txt")
pattern <- dataset[[1]]
d <- as.integer(dataset[[2]])
answer <- Neighbors(pattern, d)
output <- paste(answer, collapse = "\n")
writeClipboard(output)

cholerae <- readLines("vibrio_cholerae.txt")
