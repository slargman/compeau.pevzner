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

# 3C
pattern <- readLines("../data/rosalind_ba3c.txt")
graph <- OverlapGraph(pattern)
output <- PrintEdges(graph)
writeClipboard(output)

# 3D
dataset <- readLines("../data/rosalind_ba3d.txt")
k <- as.numeric(dataset[1])
text <- dataset[2]
graph <- DeBruijnGraph(text, k)
output <- PrintEdges(graph)
writeClipboard(output)

# 3E
patterns <- readLines("../data/rosalind_ba3e.txt")
graph <- DeBruijnGraph(patterns)
output <- PrintEdges(graph)
writeClipboard(output)

# 3F
dataset <- readLines("../data/rosalind_ba3f.txt")
graph <- AdjacencyListToGraph(dataset)
cycle <- EulerianCycle(graph)
output <- PrintPath(cycle)
writeClipboard(output)

# 3G
dataset <- readLines("../data/rosalind_ba3g.txt")
graph <- AdjacencyListToGraph(dataset)
path <- EulerianPath(graph)
output <- PrintPath(path)
writeClipboard(output)

# 3H
dataset <- readLines("../data/rosalind_ba3h.txt")
pattern <- dataset[-1]
text <- StringFromComposition(pattern)
writeClipboard(text)

# 3I
k <- as.integer(readLines("../data/rosalind_ba3i.txt"))
k_universal <- UniversalCircularString(k)
writeClipboard(k_universal)

# 3J
dataset <- readLines("../data/rosalind_ba3j.txt")
numbers <- strsplit(dataset[1], " ", fixed = T)[[1]]
d <- as.integer(numbers[2])
pattern <- dataset[-1]
text <- StringFromComposition(pattern, d = d)
writeClipboard(text)

# 3K
pattern <- readLines("../data/rosalind_ba3k.txt")
contigs <- GenerateContigs(pattern)
writeClipboard(contigs)

# 3M
adj_list <- readLines("../data/rosalind_ba3m.txt")
graph <- AdjacencyListToGraph(adj_list)
paths <- MaximalNonBranchingPaths(graph)
output <- sapply(paths, function(x) PrintPath(x, space = T))
writeClipboard(output)
