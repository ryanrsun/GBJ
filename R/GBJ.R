#' GBJ.R
#'
#' Given a vector of individual test statistics and their pairwise correlations, calculate
#' the Generalized Berk-Jones (GBJ) second-level test statistic and it's p-value.
#'
#' @param test_stats Vector of all individual (first-level) test statistics.
#' @param cor_mat A d*d matrix of the correlations between the test statistics, where
#' d is the total number of test statistics in the set.
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics.  You only need to specify EITHER cor_mat OR pairwise_cors.
#'
#' @return A list with the elements:
#' \item{GBJ}{The observed Generalized Higher Criticism test statistic.}
#' \item{GBJ_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'
#' @import stats BH
#' @importFrom Rcpp evalCpp
#' @useDynLib GBJ
#'
#' @export
#' @examples
#' # Should return statistic = 0.9248399 and p_value = 0.2670707
#' set.seed(100)
#' Z_vec <- rnorm(5) + rep(1,5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' GBJ(test_stats=Z_vec, cor_mat=cor_Z)


GBJ <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL)
{
  # Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

	# Sometimes test stats are too big for R's precision
	too_big <- which(t_vec > 8.2)
	if (length(too_big) > 0) {t_vec[too_big] <- 8.2}

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

	if (is.na(GBJ_corp) & gbj >= 20) {GBJ_corp="<1*10^(-12)"}

	return ( list(GBJ=gbj, GBJ_pvalue=GBJ_corp) )
}

