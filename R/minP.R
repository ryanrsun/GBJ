#' minP.R
#'
#' Given a vector of individual test statistics and their pairwise correlations, calculate 
#' the MinimumP (see Conneely and Boehnke, 2007) second-level test statistic and it's p-value.
#'
#' @param test_stats Vector of all individual (first-level) test statistics
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test 
#' statistics, where d is total number of test statistics in the set.  
#'
#' @return A list with the elements:
#' \item{minP}{The observed MinimumP test statistic.}
#' \item{minP_pvalue}{The p-value of this observed value, given the size of the set and 
#' correlation structure.}
#'
#' @export
#' @examples 
#' # Should return statistic = 0.05918928 and p_value = 0.2525972.
#' set.seed(100)
#' test_stats <- rnorm(5) + rep(1,5)
#' minP(test_stats, pairwise_cors=rep(0.3,10))


minP <- function(test_stats, pairwise_cors) 
{
	# Calculate minP statistic
	t_vec <- sort(abs(test_stats), decreasing=TRUE)
	d <- length(t_vec)
	minP_stat <- 1-pchisq(t_vec[1]^2, df=1)
	
	# Correct number of pairwise correlations?
	if (length(pairwise_cors) != d*(d-1)/2) {
		stop("Your pairwise correlation vector is of the wrong length!")
	}
			
	# minP bounds	 
	minP_p_bounds <- rep(minP_stat, d)
	minP_z_bounds <- qnorm(1-minP_p_bounds/2)
	minP_z_bounds <- sort(minP_z_bounds, decreasing=F)
	
	minP_z_bounds[which(minP_z_bounds > 8.2)]= 8.2
	
	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		minP_corp <- ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=rep(-999,2))	
	} else {
		minP_corp <- ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=pairwise_cors)
	}
				
	return ( list(minP=minP_stat, minP_pvalue=minP_corp) )
}

