# 5A
dataset <- readLines("../data/rosalind_ba5a.txt")
money <- as.numeric(dataset[1])
coins <- as.numeric(strsplit(dataset[2], ",", fixed = T)[[1]])
num <- DPChange(money, coins)
writeClipboard(as.character(num))

# 5B
dataset <- readLines("../data/rosalind_ba5b.txt")
nm <- as.numeric(strsplit(dataset[1], " ", fixed = T)[[1]])
n <- nm[1]
m <- nm[2]
divider <- which(dataset == "-")
down <- ReadMatrix(dataset[2:(divider - 1)])
right <- ReadMatrix(dataset[-(1:divider)])
len <- ManhattanTourist(n+1, m+1, down, right)
writeClipboard(as.character(len))
