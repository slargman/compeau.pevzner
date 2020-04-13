#' Find the minimum number of coins needed to make change
#' 
#' \code{DPChange} uses a dynamic programming approach to calculate the minimum number of coins (of denominations listed in \code{coins}) required to make change for \code{money}.
#' 
#' @param money An integer specifying the amount to be made into change.
#' @param coins A numeric vector of positive integers specifying the denominations of the coins.
#' @return The minimum number of coins with denominations \code{coins} that changes \code{money}
#' @examples
#' money <- 40
#' coins <- c(1, 5, 10, 20, 25, 50)
#' DPChange(money, coins)
DPChange <- function(money, coins){
	min_num_coins <- rep(Inf, money + 1)
	min_num_coins[1] <- 0
	for (m in 2:(money + 1)) {
		for (i in 1:length(coins)) {
			if (m - 1 >= coins[i]) {
				if (min_num_coins[m - coins[i]] + 1 < min_num_coins[m]) {
					min_num_coins[m] <- min_num_coins[m - coins[i]] + 1
				}
			}
		}
	}
	return(min_num_coins[money + 1])
}

rows <- c("1 0 2 4 3", "4 6 5 2 1", "4 4 5 2 1", "5 6 8 5 3") 
ReadMatrix(rows)
ReadMatrix <- function(rows){
	mat <- strsplit(rows, split = " ", fixed = T)
	mat <- t(sapply(mat, as.numeric))
	return(mat)
}

n <- 4
m <- 4
rows1 <- c("1 0 2 4 3", "4 6 5 2 1", "4 4 5 2 1", "5 6 8 5 3") 
down <- ReadMatrix(rows1)
rows2 <- c("3 2 4 0", "3 2 4 2", "0 7 3 3", "3 3 0 2", "1 3 2 2")
right <- ReadMatrix(rows2)
ManhattanTourist(n+1, m+1, down, right)
ManhattanTourist <- function(n, m, down, right){
	s <- matrix(nrow = n, ncol = m)
	s[1, 1] <- 0
	# fill in left edge
	for (i in 2:n) {
		s[i, 1] <- s[i - 1, 1] + down[i - 1, 1]
	}
	# fill in top edge
	for (j in 2:m) {
		s[1, j] <- s[1, j - 1] + right[1, j - 1]
	}
	# fill in the rest
	for (i in 2:n) {
		for (j in 2:m) {
			s[i, j] <- max(s[i - 1, j] + down[i - 1, j], s[i, j - 1] + right[i, j - 1])
		}
	}
	return(s[n,m])
}

v <- "AACCTTGG"
w <- "ACACTGTGA"
LCSBacktrack(v, w)
LCSBacktrack <- function(v,w){
	length_v <- nchar(v)
	length_w <- nchar(w)

	vec_v <- strsplit(v, split = "")[[1]]
	vec_w <- strsplit(w, split = "")[[1]]

	s <- matrix(nrow = length_v + 1, ncol = length_w + 1)
	backtrack <- matrix(nrow = length_v + 1, ncol = length_w + 1)
	
	# fill in left edge
	for (i in 1:(length_v + 1)){
		s[i, 1] <- 0
		backtrack[i, 1] <- "v"
	}

	# fill in top edge
	for (j in 1:(length_w + 1)){
		s[1, j] <- 0
		backtrack[1, j] <- ">"
	}

	# fill in the rest
	for (i in 2:(length_v + 1)){
		for (j in 2:(length_w + 1)){

			# fill in s[i,j]
			s[i, j] <- max(s[i - 1, j], s[i, j - 1], s[i - 1, j - 1] + (vec_v[i - 1] == vec_w[j - 1])) 

			# fill in backtrack[i,j]
			if (s[i, j] == s[i - 1, j]) {
				backtrack[i, j] <- "v"
			} else if (s[i, j] == s[i, j - 1]) {
				backtrack[i, j] <- ">"
			} else if ((s[i,j] == (s[i - 1, j - 1] + 1)) && (vec_v[i - 1] == vec_w[j - 1])){
				backtrack[i, j] <- "L"
			}
		}
	}

	return(backtrack)
}

v <- "AACCTTGG"
w <- "ACACTGTGA"
backtrack <- LCSBacktrack(v, w)
OutputLCS(backtrack, v, nchar(v), nchar(w))
OutputLCS <- function(backtrack, v, i, j){
	if ((i == 0) || (j == 0)) {
		return("")
	}

	if (backtrack[i + 1, j + 1] == "v"){
		OutputLCS(backtrack, v, i - 1, j)
	} else if (backtrack[i + 1, j + 1] == ">"){
		OutputLCS(backtrack, v, i, j - 1)
	} else {
		return(paste0(OutputLCS(backtrack, v, i - 1, j - 1), substring(v, i, i)))
	}
}
