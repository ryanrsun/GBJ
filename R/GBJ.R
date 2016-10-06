#' GBJ.R
#'
#' Given a vector of individual test statistics and their pairwise correlations, calculate 
#' the Generalized Berk-Jones (GBJ) second-level test statistic and it's p-value.
#'
#' @param test_stats Vector of all individual (first-level) test statistics
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test 
#' statistics, where d is total number of test statistics in the set.  
#'
#' @return A list with the elements:
#' \item{GBJ}{The observed Generalized Higher Criticism test statistic.}
#' \item{GBJ_pvalue}{The p-value of this observed value, given the size of the set and 
#' correlation structure.}
#'
#' @import stats
#' @importFrom Rcpp evalCpp
#' @useDynLib GBJ
#'
#' @export
#' @examples 
#' # Should return statistic = 0.9248399 and p_value = 0.2670707
#' set.seed(100)
#' test_stats <- rnorm(5) + rep(1,5)
#' GBJ(test_stats, pairwise_cors=rep(0.2,10))


GBJ <- function(test_stats, pairwise_cors) 
{
	# Ensure that the thresholds are sorted in descending order, largest first.
	t_vec <- sort(abs(test_stats), decreasing=TRUE)
	d <- length(t_vec)
	
	# Correct number of pairwise correlations?
	if (length(pairwise_cors) != d*(d-1)/2) {
		stop("Your pairwise correlation vector is of the wrong length!")
	}

	# Calculate the observed GBJ statistic
	GBJ_stats <- GBJ_objective(t_vec=t_vec, d=d, pairwise_cors=pairwise_cors)
	gbj <- max(GBJ_stats)
	
	# Calculate p_value
	GBJ_p_list <- GBJ_pvalue(observed_gbj=gbj, d=d, pairwise_cors=pairwise_cors)
	GBJ_corp=GBJ_p_list$GBJ_corp
				
	return ( list(GBJ=gbj, GBJ_pvalue=GBJ_corp) )
}

