# 2A
dataset <- readLines("../data/rosalind_ba2a.txt")
k <- as.integer(substring(dataset[1], 1, 1)[[1]])
d <- as.integer(substring(dataset[1], 3, 3)[[1]])
dna <- dataset[-1]
answer <- MotifEnumeration(dna, k, d)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 2B
dataset <- readLines("../data/rosalind_ba2b.txt")
k <- as.integer(dataset[1])
dna <- dataset[-1]
answer <- MedianString(dna, k)
output <- answer
writeClipboard(output)

# 2C
dataset <- readLines("../data/rosalind_ba2c.txt")
text <- dataset[1]
k <- as.integer(dataset[2])
A <- as.numeric(strsplit(dataset[3], split = " ")[[1]])
C <- as.numeric(strsplit(dataset[4], split = " ")[[1]])
G <- as.numeric(strsplit(dataset[5], split = " ")[[1]])
T <- as.numeric(strsplit(dataset[6], split = " ")[[1]])
profile <- rbind(A, C, G, T)
answer <- ProfileMostProbableString(text, k, profile)
output <- paste(answer, collapse = " ")
writeClipboard(output)

# 2D
dataset <- readLines("../data/rosalind_ba2d.txt")
numbers <- as.numeric(strsplit(dataset[1], split = " ")[[1]])
k <- numbers[1]
t <- numbers[2]
dna <- dataset[-1]
answer <- GreedyMotifSearch(dna, k, t)
output <- MotifString(answer)
writeClipboard(output)

#2E
dataset <- readLines("../data/rosalind_ba2e.txt")
numbers <- as.numeric(strsplit(dataset[1], split = " ")[[1]])
k <- numbers[1]
t <- numbers[2]
dna <- dataset[-1]
answer <- GreedyMotifSearchPseudocounts(dna, k, t)
output <- MotifString(answer)
writeClipboard(output)

#2F
dataset <- readLines("../data/rosalind_ba2f.txt")
numbers <- as.numeric(strsplit(dataset[1], split = " ")[[1]])
k <- numbers[1]
t <- numbers[2]
dna <- dataset[-1]
answer <- RandomizedMotifSearch(dna, k, t, iter = 300)
output <- MotifString(answer)
writeClipboard(output)

#2G
dataset <- readLines("../data/rosalind_ba2g.txt")
numbers <- as.numeric(strsplit(dataset[1], split = " ")[[1]])
k <- numbers[1]
t <- numbers[2]
N <- numbers[3]
dna <- dataset[-1]
answer <- GibbsSampler(dna, k, t, N, iter = 9)
output <- MotifString(answer)
writeClipboard(output)
