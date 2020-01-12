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
