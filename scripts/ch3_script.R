# 3A
dataset <- readLines("../data/rosalind_ba3a.txt")
k <- as.numeric(dataset[1])
text <- dataset[2]
answer <- StringComposition(text, k)
output <- answer
writeClipboard(output)

# 3B
pattern <- readLines("../data/rosalind_ba3b.txt")
answer <- StringFromGenomePath(pattern)
writeClipboard(answer)
