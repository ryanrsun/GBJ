// ebb_crossprob_cor.c
// Ryan Sun, Harvard University
// July 17, 2016

// This binary calculates an approximation to the p-value of GOF tests.

// Usage is ./ebb_crossprob_cor [NUM_STATS] [BOUNDS_FILE] [CORRELATIONS_FILE]
// where NUM_STATS is the size of your set, BOUNDS_FILE is the name of the file
// holding the boundaries (0<=b_1<=...<=b_NUM_STATS) which come from inversion
// of the GOF statistic, and CORRELATIONS_FILE is the name of the file holding
// all the NUM_STATS*(NUM_STATS-1)/2 pairwise correlations between the observations
// in your set.
// Both BOUNDS_FILE and CORRELATIONS_FILE should be a single column of numbers
// with no headers/labels/blank spaces/other characters.

// One difference from before is that we calculate the conditional moments upfront in
// a separate routine.
// We are also now using the EBB distribution of Prentice (1988) instead of the standard
// Beta-Binomial.

// The routines below are:
// (1) avg_cond_covar - Calculate the average conditional covariance between any two indicators
//  Y_i,Y_j where Y_i=P(|Z_i|>=t_k).
// (2) match_moments - Given an average conditional covariance and a conditional mean, calculate
// \lambda and \gamma for the EBB PMF.  Also check to ensure both are inside the allowable parameter space.
// (3) eval_EBB_PMF - Self-explanatory
// (4) calc_qka - Called to fill each entry of the d*d matrix leading to final
// p-value.  Sums over all m>=a in P(S(t_k)=a|S(t_k)=m).
// (5) calc_allq - Loop through each entry of the d*d matrix until we get the (d,1) element.

// Need this for the boost functions on windows machines
// [[Rcpp::depends(BH)]]
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <numeric>
#include <boost/math/special_functions/erf.hpp>		// For normal PDF/CDF
#include <boost/math/constants/constants.hpp>			// For pi()
#include <Rcpp.h>

using namespace Rcpp;


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


// To square a number faster than pow().
inline double fast_square(const double & x)
{
    return (x*x);
}

// To cube a number faster than pow().
inline double fast_cube(const double & x)
{
    return (x*x*x);
}

// PDF of a standard normal RV.
inline double dnorm(const double & x)
{
    return (exp(-fast_square(x) / 2.0) /
            sqrt(2.0*boost::math::constants::pi<long double>()));
}

// Survival function of a standard normal, 1-F_x(x)
inline double surv(const double & x)
{
    return (1 - 0.5 * erfc(-x/sqrt(2.0)));
}



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// For the independence case, good to have lookup table of log-factorials
// Assume d doesn't go past the limit of 'int,' which is 32,767.
bool create_logftable(const int & d, std::vector<double> & log_ftable)
{
    // 0! = 1! = 1
    log_ftable[0] = 0.0;
    log_ftable[1] = 0.0;

    for (int iii=2; iii<=d; ++iii)
    {
        log_ftable[iii] = lgamma(iii+1);
    }

    return 0;
}

// Create a d(d-1)/2 by 5 table of each \rho_ij to the 2/4/6/8/10 power.
// These are reused d times in the calculation of the variance of S(t).
bool create_rtable(const int & d,
                   const std::vector <double> & r_vec,
                   std::vector< std::vector<double> > & r_table)
{
    long double num_cor = d*(d-1)/2;
    std::vector <double> r_two(num_cor);
    std::vector <double> r_four(num_cor);
    std::vector <double> r_six(num_cor);
    std::vector <double> r_eight(num_cor);
    std::vector <double> r_ten(num_cor);

    // Long is approximately 2*10^9 which is enough to cover up to d~40000.
    for (long iii=0; iii<num_cor; ++iii)
    {
        r_two[iii] = fast_square(r_vec[iii]);
        r_four[iii] = pow(r_vec[iii], 4.0);
        r_six[iii] = pow(r_vec[iii], 6.0);
        r_eight[iii] = pow(r_vec[iii], 8.0);
        r_ten[iii] = pow(r_vec[iii], 10.0);
    }

    r_table[0] = r_two;
    r_table[1] = r_four;
    r_table[2] = r_six;
    r_table[3] = r_eight;
    r_table[4] = r_ten;

    return 0;
}

// Create a d*5 matrix holding the values of the 2/4/6/8/10th Hermite polynomials
// evaluated at t_k from k=1,...,d.
// Each value will be reused d(d-1)/2 times in the calculation of the variance of S(t).
// We also pre-divide each of the He_(2n-1) terms by (2n!) to save a little more time.
bool create_hermtable(const int & d,
                      const std::vector <double> & t_vec,
                      std::vector<std::vector<double> > & herm_table)
{
    // Hold each of the five Hermite polynomial terms.
    std::vector <double> herm_one(d);
    std::vector <double> herm_two(d);
    std::vector <double> herm_three(d);
    std::vector <double> herm_four(d);
    std::vector <double> herm_five(d);

    for (int iii=0; iii<d; ++iii)
    {
        herm_one[iii] = fast_square( t_vec[iii] ) / 2.0;
        herm_two[iii] = fast_square( fast_cube(t_vec[iii]) - 3.0*t_vec[iii] ) / 24.0;
        herm_three[iii] = fast_square( pow(t_vec[iii], 5) - 10.0*fast_cube(t_vec[iii]) +15.0*t_vec[iii] ) / 720.0;
        herm_four[iii] = fast_square( pow(t_vec[iii], 7) - 21.0*pow(t_vec[iii], 5)
                                     + 105.0*fast_cube(t_vec[iii]) - 105.0*t_vec[iii] ) / 40320.0;
        herm_five[iii] = fast_square( pow(t_vec[iii], 9) - 36.0*pow(t_vec[iii], 7)
                                     + 378.0*pow(t_vec[iii], 5) - 1260.0*fast_cube(t_vec[iii]) + 945.0*t_vec[iii] ) / 3628800.0;
    }

    // Populate first dimension.
    herm_table[0] = herm_one;
    herm_table[1] = herm_two;
    herm_table[2] = herm_three;
    herm_table[3] = herm_four;
    herm_table[4] = herm_five;

    return 0;
}


// Calculate the average conditional covariance term, P(|Z_i|,|Z_j| >= t_k)
// for k=1,...,d.
bool avg_cond_covar(const int & d,
                    const std::vector<double> & t_vec,
                    const std::vector< std::vector<double> > & r_table,
                    const std::vector< std::vector<double> > & herm_table,
                    std::vector< double > & covar_vec)
{
    // Numerator and denominator vectors for the tricky infinite sum in the variance.
    double num_cor = d*(d-1)/2;
    double cond_mean_sq = 0.0;
    std::vector<double> numerator_vec(num_cor);
    std::vector<double> denominator_vec(num_cor);

    // First the unconditional covariance at k=1.
    double sum_term = 0.0;
    double surv_term = 4 * fast_square( surv(t_vec[0]) );
    double phi_term = 4 * fast_square( dnorm(t_vec[0]) );
    for (long cor_it=0; cor_it<num_cor; ++cor_it)
    {
        numerator_vec[cor_it] = surv_term + phi_term *
        (r_table[0][cor_it]*herm_table[0][0] +
         r_table[1][cor_it]*herm_table[1][0] +
         r_table[2][cor_it]*herm_table[2][0] +
         r_table[3][cor_it]*herm_table[3][0] +
         r_table[4][cor_it]*herm_table[4][0]);
        sum_term += numerator_vec[cor_it];
    }
    cond_mean_sq = fast_square( 2*surv(t_vec[0]) );
    covar_vec[0] = ( sum_term - num_cor*cond_mean_sq ) / num_cor;

    // Now calculate the conditional p_n(1) and p_n(2) at k=2,...,d.
    for (int kkk=2; kkk<=d; ++kkk)
    {
        // For a test like GBJ the smallest half t_k are the same, so we don't need to call
        // calc_qka or calculate conditional moments for first half of rows.
        // Make sure our tolerance here 10^-8 is the same as in calc_allq.
        if ( fabs(t_vec[kkk-1] - t_vec[kkk-2]) < pow(10.0, -8.0) )
        {
            covar_vec[kkk-1] = covar_vec[kkk-2];
            continue;
        }

        // The conditional variance at t_k.
        denominator_vec.swap(numerator_vec);
        sum_term = 0.0;
        surv_term = 4 * fast_square( surv(t_vec[kkk-1]) );
        phi_term = 4 * fast_square( dnorm(t_vec[kkk-1]) );
        for (long cor_it=0; cor_it<num_cor; ++cor_it)
        {
            numerator_vec[cor_it] = surv_term + phi_term *
            (r_table[0][cor_it]*herm_table[0][kkk-1] +
             r_table[1][cor_it]*herm_table[1][kkk-1] +
             r_table[2][cor_it]*herm_table[2][kkk-1] +
             r_table[3][cor_it]*herm_table[3][kkk-1] +
             r_table[4][cor_it]*herm_table[4][kkk-1]);
            sum_term += numerator_vec[cor_it] / denominator_vec[cor_it];
        }

        // Record
        cond_mean_sq = fast_square( surv(t_vec[kkk-1]) / surv(t_vec[kkk-2]) );
        covar_vec[kkk-1] = ( sum_term - num_cor*cond_mean_sq ) / num_cor;
    }

    return 0;
}

// Evaluate the EBB PMF for a range of n, so we don't have to repeat
// multiplications for Pr[S(t_k)=a|S(t_k-1)=m] for m=a:(d-k+1).  Here 'y' is a.
bool eval_EBB_PMF_allN(const int & max_n,
                       const int & y,                  // min_n = y
                       const double & lambda,
                       const double & gamma,
                       std::vector<double> & PMF_vec)
{
    // If (d-k+1) <=1 then we don't need to do these calculations.
    if (max_n < 2) { return 0;}

    double prob_mass = 1.0;
    double min_n;
    if (y < 2)        // Here y=0/1 and max_n >= 2; can fill the 0/1 slots of PMF with nonsense.
    {
        PMF_vec[0] = 0.0;
        PMF_vec[1] = 0.0;

        if (y == 0)
        {
            prob_mass = (1-lambda) * (1-lambda+gamma) / (1+gamma);      // For a=0 and m=2:(d-k+1)
        } else if (y == 1)
        {
            prob_mass = lambda * (1-lambda) / (1+gamma);                // For a=1 and m=2:(d-k+1)
        }

        // Start the next phase at y=2.
        min_n = 2;
        PMF_vec[2] = prob_mass;
    }
    else
    {
        // Here we just do one calculation and finish.
        for (int iii=0; iii<y; ++iii)
        {
            prob_mass = prob_mass * (lambda+gamma*iii) / (1+gamma*iii);
        }
        PMF_vec[y] = prob_mass;
        min_n = y;
    }

    if (min_n == max_n)             // Done if a=(d-k+1).
    {
        return 0;
    }

    // Not done, a<(d-k+1).
    int n_diff = max_n - min_n;
    for (int jjj=1; jjj<=n_diff; ++jjj)             // One multiplication and store for each j.
    {
        prob_mass = prob_mass * (1-lambda+gamma*(min_n-y-1+jjj)) / (1+gamma*(min_n-1+jjj));
        PMF_vec[min_n+jjj] = prob_mass;
    }

    return 0;
}

// Calculate q_k,a by summing over q_k,a|S(t)=m for m=a:(d-k)
double calc_qka(const int & d,
                const int & k,
                const int & a,
                const std::vector<double> & prev_row,
                const std::vector<double> & log_ftable,
                const bool & ind_flag,
                const double & lambda,
                const double & gamma)
{
    double q_ka = 0.0;
    double min_gamma;
    double m_choose_a;
    std::vector<double> PMF_vec(d+1);

    // To hold the the EBB probability (w/o factorial part) for m=a:(d-k+1)
    if (!ind_flag)
    {
        bool all_EBB_status = eval_EBB_PMF_allN((d-k+1),
                                                a,
                                                lambda,
                                                gamma,
                                                PMF_vec);

         if (all_EBB_status)
    	{
    		return(-1);
   	 	}
    }

    // Sum over all possible values of S(t)=m
    for (int mmm=a; mmm<=(d-k+1); ++mmm)
    {
        m_choose_a = std::exp( log_ftable[mmm] - log_ftable[a] - log_ftable[mmm-a] );

        // Use binomial if m=0/1 or if independence flag or if gamma outside parameter space.
        min_gamma = std::max(-lambda/(mmm-1), -(1-lambda)/(mmm-1));
        if (mmm<=1 || ind_flag || gamma<min_gamma)
        {
            q_ka += prev_row[mmm] * m_choose_a * pow(lambda, a) * pow((1.0-lambda), (mmm-a));
        }
        else          // EBB PMF
        {
            q_ka += prev_row[mmm] * m_choose_a * PMF_vec[mmm];
        }
    }

    return q_ka;
}

// The p-value calculation 'master' function, the interface between the math and main fn.
// Loop through all k, calculatin q_k,a from a=0:(t-k) until we get to q_d,0/
// Return the p-value.
double calc_allq(const int & d,
                 const std::vector <double> & t_vec,
                 const std::vector <double> & r_vec,
                 const bool & ind_flag)
{
    // Fill the d(d-1)/2 by 5 table of correlations to the 2/4/6/8/10 power.
    // Fill the d*5 table of Hermite polynomial terms in the S(t) conditional variance calculation.
    long num_cor = d*(d-1)/2;
    std::vector<std::vector<double> > r_table(5, std::vector<double>(num_cor));
    std::vector<std::vector<double> > herm_table(5, std::vector<double>(d));
    bool filled_hermtable = create_hermtable(d,
                                             t_vec,
                                             herm_table);
    bool filled_rtable = create_rtable(d,
                                       r_vec,
                                       r_table);

    if (filled_hermtable || filled_rtable)
    {
        std::cout << "Problem creating lookup tables." << std::endl;
        exit(1);
    }

    // Fill the d*1 vector of log-factorials
    std::vector<double> log_ftable(d+1);
    bool filled_logftable = create_logftable(d,
                                             log_ftable);
    if (filled_logftable)
    {
    	return(-1);
    }

    // Calculate the conditional average pairwise correlation for k=1,...,d.
    std::vector<double> covar_vec(d);
    bool filled_condcovar = avg_cond_covar(d,
                                           t_vec,
                                           r_table,
                                           herm_table,
                                           covar_vec);
    if (filled_condcovar)
    {
        std::cout << "Problem calculating moments." << std::endl;
        exit(1);
    }

    // We don't actually need to hold the entire d*d matrix, just have to always
    // know the previous row.
    std::vector<double> prev_row(d+1);
    std::vector<double> current_row((d+1), 0.0);

    // Initial conditions.
    double prev_bound = 0.0;
    current_row[d] = 1.0;

    // Loop through k=1,...,d.
    double lambda = 1.0;
    double gamma = 1.0;
    double rho;
    double avg_cond_covar;
    double temp_qka;
    for (int kkk=1; kkk<=d; ++kkk)
    {
        prev_row = current_row;

        // If new bound same as previous bound (make sure tolerance is same as match_moments).
        if (fabs(t_vec[kkk-1] - prev_bound) < pow(10.0, -8.0))
        {
            current_row[d-kkk+1] = 0.0;
            prev_bound = t_vec[kkk-1];
            continue;
        }

        // New bound, new probabilities for the row.
        std::fill(current_row.begin(), current_row.end(), 0.0);

        // Match moments once for each row
        lambda = surv(t_vec[kkk-1]) / surv(prev_bound);
        avg_cond_covar = covar_vec[kkk-1];
        rho = avg_cond_covar / (lambda*(1-lambda));
        gamma = rho / (1-rho);

        // For each k, we want q_k,a for a=0:(t-k).
        for (int aaa=0; aaa<=(d-kkk); ++aaa)
        {
            temp_qka = calc_qka(d,
                                kkk,
                                aaa,
                                prev_row,
                                log_ftable,
                                ind_flag,
                                lambda,
                                gamma);
            current_row[aaa] = temp_qka;

        }

        // Update so we can check if same as next bound, also for matching lambda.
        prev_bound = t_vec[kkk-1];
    }

    return (1-current_row[0]);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Main function reads in our boundary points and correlations, the
// runs calc_allq and prints the p-value.
// If the `correlations` argument is less than -1, then this will trigger
// the independence flag.

// [[Rcpp::export]]

double ebb_crossprob_cor_R(int d, NumericVector bounds, NumericVector correlations) {

    // For the Rcpp version, just input the bounds as NumericVectors and transfer them
    // to std::vectors.
    long num_cor = d*(d-1)/2;
    std::vector<double> boundary_pts;
    std::vector<double> cors_vec;

    // Reserve memory to prevent fragmentation, easy since we know the exact size
    boundary_pts.reserve(d);

    // Put the boundary pts into array, one by one.
    // Should be sorted in increasing order.
    for (int iii=0; iii<d; ++iii)
    {
        boundary_pts.push_back(bounds[iii]);
    }

    // If the last argument is -999 then no need for correlation vector
    // Carry through the independence flag so that we don't ever need to matchMoments()
    bool indFlag = false;
    if (correlations[1] < -1)
    {
        indFlag = true;
        cors_vec.push_back(-1.0);
    }
    else
    {
        // Put the correlations into array, one by one.
        cors_vec.reserve(num_cor);
        for (long iii=0; iii<num_cor; ++iii)
        {
            cors_vec.push_back(correlations[iii]);
        }
    }

    // Calculate the p-value.
    double pvalue = calc_allq(d, boundary_pts, cors_vec, indFlag);

    // Return the p-value.
    return pvalue;
}




