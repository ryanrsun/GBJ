#' HC.R
#'
#' Given a vector of individual test statistics and their pairwise correlations, calculate
#' the Higher Criticism second-level test statistic and it's p-value.
#'
#' @param test_stats A scalar or vector of threshold points (magnitudes of test statistics)
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics, where d is total number of test statistics in the set.
#'
#' @return A list with the elements:
#' \item{HC}{The observed Higher Criticism test statistic.}
#' \item{HC_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'
#' @export
#' @examples
#' # Should return statistic = 2.067475 and p_value = 0.2755146
#' set.seed(100)
#' test_stats = rnorm(5) + rep(1,5)
#' HC(test_stats, pairwise_cors=rep(0.2,10))

HC <- function(test_stats, pairwise_cors) {

	# Ensure that the thresholds are sorted in descending order, largest first.
	t_vec <- sort(abs(test_stats), decreasing=TRUE)
	d <- length(t_vec)

	# Sometimes test stats are too big for R's precision
	too_big <- which(t_vec > 8.2)
	if (length(too_big) > 0) {t_vec[too_big] <- 8.2}

	# Correct number of pairwise correlations?
	if (length(pairwise_cors) != d*(d-1)/2) {
		stop("Your pairwise correlation vector is of the wrong length!")
	}

	# Calculate HC objectives
	p_values <- 1-pchisq(t_vec^2, df=1)
	i_vec <- 1:d
	HC_stats <- sqrt(d) * (i_vec/d - p_values) / sqrt(p_values*(1-p_values))

	# Observed HC statistic
	h <- max(HC_stats)

	# Calculate p-value
	if (h<=0) {
		return ( list(HC=0, HC_pvalue=1) )
	}

	# BJ bounds
	HC_p_bounds <- rep(NA, d)

	# Explicit inverse of HC to find the p-value bounds
	HC_p_bounds <- ((2*i_vec+h^2)/d - sqrt((2*i_vec/d+h^2/d)^2 - 4*i_vec^2/d^2 - 4*i_vec^2*h^2/d^3))/(2*(1+h^2/d))
	HC_z_bounds <- qnorm(1-HC_p_bounds/2)
	HC_z_bounds <- sort(HC_z_bounds, decreasing=F)

	# qnorm can't handle more precision than 10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	HC_z_bounds[which(HC_z_bounds > 8.2)]= 8.2

	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		HC_corp <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=rep(-999,2))
	} else {
		HC_corp <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=pairwise_cors)
	}


	return ( list(HC=h, HC_pvalue=HC_corp) )
}
