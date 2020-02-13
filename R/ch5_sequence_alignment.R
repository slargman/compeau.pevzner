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
