#' calc_score_stats.R
#'
#' Starting with individual-level data on p factors, generate score test statistics for each
#' factor for input into GBJ/GHC/HC/BJ/minP.  Also get the correlations between these test statistics. 
#' Designed to be used with linear or logistic regression null models.
#'
#' @param null_model A fitted R regression model
#' @param factor_matrix An n*p matrix with each factor as one column.  There should be no missing data.
#' @param model_type Either "linear" or "logistic"
#'
#' @return A list with the elements:
#' \item{test_stats}{The p score test statistics.}
#' \item{cor_mat}{The p*p matrix giving the pairwise correlation of every two test statistics.}
#'
#' @export
#' @examples 
#' Y <- rbinom(n=100, size=1, prob=0.5)
#' null_mod <- glm(Y~1, family=binomial(link="logit"))
#' factor_mat <- matrix(data=rnorm(n=100*5), nrow=100)
#' calc_score_stats(null_mod, factor_mat, "logistic")

calc_score_stats <- function(null_model, factor_matrix, model_type) {
	
	X_mat <- model.matrix(null_model)
	fitted_Y <- null_model$fitted.values
	actual_Y <- null_model$y
	
	# Only difference between linear and logistic procedure
	if (model_type == 'logistic') {
		W_vec <- fitted_Y * (1-fitted_Y)
	} else if (model_type == 'linear') {
		W_vec <- rep(summary(null_model)$sigma^2, nrow(X_mat))
	} else {
		stop("Invalid model type")
	}
	
	########################
	# EZ Mode if linear regression, no additional covariates except for intercept
	if (model_type == 'linear' & formula(null_model) == "Y ~ 1") {
		num_sub <- nrow(X_mat)
		sig_sq_hat <- sum( (actual_Y - fitted_Y)^2 ) / (num_sub-1)
		test_stats <- rep(NA, d) 
		denominators <- rep(NA, d)
		for(kkk in 1:d) 
		{
			tempF<- factor_matrix[,kkk]
			score_num <- t(tempF) %*% (actual_Y-fitted_Y) 
			score_denom <- sqrt(sig_sq_hat * (sum(tempF^2) - mean(tempF)^2*num_sub))
			denominators[kkk] <- score_denom
			test_stats[kkk] <- score_num / score_denom
		}
		est_cor <- cor(factor_matrix)
		
		# Return from here
		return ( list(test_stats=test_stats, cor_mat=est_cor) )	
	}
	
	########################
	# Regular mode
	W_mat <- diag(W_vec)
	P_mat <- W_mat - W_mat%*%X_mat %*% solve(t(X_mat)%*%W_mat%*%X_mat) %*% t(X_mat)%*%W_mat
	
	# Now our score test
	d <- ncol(factor_matrix)
	test_stats <- rep(NA, d)
	denominators <- rep(NA, d)
	for(kkk in 1:d) 
	{
		# Pick out next SNP, conduct score test (no additional covariates).
		tempF <- factor_matrix[,kkk]
		score_num <- t(tempF) %*% (actual_Y-fitted_Y)
		score_denom <- sqrt(tempF %*% P_mat %*% tempF)
		test_stats[kkk] <- score_num / score_denom
		denominators[kkk] <- score_denom
	}

	# Estimate the correlation matrix for the test statistics.
	# Same as the correlation matrix of the SNPs if linear regression and no additional covariates.
	est_cor <- matrix(data=NA, nrow=d, ncol=d)
	for (temp_row in 2:d)
	{
		for (temp_col in 1:(temp_row-1))
		{
			est_cor[temp_row, temp_col] <- t(factor_matrix[,temp_row]) %*% P_mat %*% factor_matrix[,temp_col] / 										(denominators[temp_row] * denominators[temp_col])
			est_cor[temp_col, temp_row] <- est_cor[temp_row, temp_col]
		}
	}
	
	return ( list(test_stats=test_stats, cor_mat=est_cor) )	
}
