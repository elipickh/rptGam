#' @title Estimate the number of BCa bootstrap replicates 
#'
#' @param A A numeric value, indicating the acceleration factor. Can be estimated from the data 
#' (i.e., using jackknife) or, if A is bnull, will be estimated from boots as in coxed::bca, 
#' but using the point estimate instead of mean of boots in the calcualtion.
#' @param point A numeric value. The sample statistic from the original dataset.
#' @param perc_dev A numeric in (0, 100), indicating the maximum percentage deviation of each 
#' of the upper and lower confidence interval lengths. The lower this value, the more boostraps 
#' are typically required.
#' @param prob: A numeric n (0,1), indicating the approximate probability that perc_dev holds.
#' @details 
#' "The accuracy obtained by a given choice of B is random, because the bootstrap simulations are random.
#'  To determine an appropriate value of B, we specify a bound on the percentage deviation, denoted pdb, 
#'  and we require that the actual percentage deviation is less than this bound with a specified probability, 
#'  1 - r, close to one. 
#'  The three-step method takes pdb and r as given and specifies a data-dependent method 
#'  of determining a value of B...such that the desired level of accuracy is achieved. For example, one might 
#'  take (pdb, r) = (10, .05). In this case, the three-step method determines a value B* such that the 
#'  percentage deviation of the upper and lower confidence interval lengths is less than 10% each with 
#'  approximate probability .95." (Andrews & Buchinsky, 2002)
#' Authors recommend iterating through this process for (possibly) better finite accuracy.
#' @return An integer, indicating the estimated number of bootstrap replicates. This value is either the 
#' initial estimate (if boots is null) or the final estimate (if boots is not null). 
#' 
#' @examples
#'
#' @references Andrews, D. W., & Buchinsky, M. (2002). On the number of bootstrap repetitions 
#' for BCa confidence intervals. Econometric Theory, 18(4), 962-984.
#'
#' @author Eliezer Pickholtz (eyp3@@cornell.edu)
#' @export
#'
bca_coverage <- function(A, ci = 0.95, perc_dev = 10, prob = 0.95, boots = NULL, point = NULL) {
  a = (1 - ci)/2
  # a should be between 0.01, 0.99, according to paper
  if (a > 0.99 | a < 0.01) stop('a should be in [0.01, 0.99]')
  r = 1 - prob
  pdb = perc_dev

  if (is.null(boots)) {
    B1 = 
        # according to formula:
        ceiling(
        # according to paper examples:
        #floor(
        10000 * 
        (a * (1 - a) - 
        2 * a * dnorm(qnorm(a)) / dnorm(0) + 
        (dnorm(qnorm(a))^2) / (dnorm(0)^2)) *
        (qnorm(1 - r/2)^2) /
        ((qnorm(a) * dnorm(qnorm(a)) * pdb)^2)
        )
    return(B1)
  } 
  B1 = length(boots)
  # initial bias correction:
  Z1 = qnorm(sum(boots < point) / B1)
  if (is.null(A)) {
    L = (B1 - 1) * (point - boots);
    A = sum(L^3) / (6 * sum(L^2)^1.5)
  }
  a1l = max(pnorm(Z1 + ((Z1 + qnorm(a)) / (1 - A * (Z1 + qnorm(a))))), 0.01)
  a1u = min(pnorm(Z1 + ((Z1 + qnorm(1 - a)) / (1 - A * (Z1 + qnorm(1 - a))))), 0.99)
  v1l = max(floor((B1 + 1) * a1l), 1)
  v1u = min(B1, ceiling((B1 + 1) * a1u))
  C_func = function(a) {((1.5 * (qnorm(1 - a/2)^2) * (dnorm(qnorm(1 - a))^2)) / 
                       (2 * (qnorm(1 - a)^2) + 1)) ^ 
                      (1/3)}
  m1l = ceiling(C_func(a1l) * (B1^(2/3)))
  m1u = ceiling(C_func(1 - a1u) * (B1^(2/3)))
  B2l = 
    ceiling(
  10000 * (a * (1 - a) - 2 * a * dnorm(qnorm(a)) / dnorm(0) + (dnorm(qnorm(a))^2) / (dnorm(0)^2)) * 
  (qnorm(1 - r / 2)^2) *
  #((B1 / (2 * m1l))^2) * ((boots[v1l + m1l] - boots[v1l - m1l])^2) / (((point - boots[v1l]) * pdb)^2)
      ((B1 / (2 * m1l))^2) * ((boots[min(B1, v1l + m1l)] - boots[max(1, v1l - m1l)])^2) / (((point - boots[v1l]) * pdb)^2)
        )
  B2u = 
      ceiling(
    10000 * (a * (1 - a) - 2 * a * dnorm(qnorm(a)) / dnorm(0) + (dnorm(qnorm(a))^2) / (dnorm(0)^2)) * 
    (qnorm(1 - r / 2)^2) *
    #((B1 / (2 * m1u))^2) * ((boots[v1u + m1u] - boots[v1u - m1u])^2) / (((boots[v1u] - point) * pdb)^2)
      ((B1 / (2 * m1u))^2) * ((boots[min(B1, v1u + m1u)] - boots[max(1, v1u - m1u)])^2) / (((boots[v1u] - point) * pdb)^2)
          )

  B = max(B1, B2l, B2u)
  B
  
}