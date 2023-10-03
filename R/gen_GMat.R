#' gen_GMat.R
#'
#' Generate the multivariate binary random genotype matrix for simulation studies.
#'
#' @param n the number of observations to generate (sample size).
#' @param q the number of SNPs.
#' @param rho the correlation for the binary correlation matrix.
#'
#' @return A n*q correlated matrix for simulated SNPs.
#'
#' @importFrom stats toeplitz
#' @importFrom bindata rmvbin
#'
#' @export
#' @examples
#' n <- 1000
#' q <- 20
#' rho <- 0.1
#' gen_GMat(n, q, rho)
#'
#'
gen_GMat <- function(n, q, rho){
  # Construct a binary correlation matrix
  cmat <- stats::toeplitz(c(1, rep(rho, q - 1))) # q * q
  meanparam1 <- runif(q, .01, .05)
  meanparam2 <- runif(q, .01, .05)
  x <- suppressWarnings(bindata::rmvbin(n, margprob = meanparam2, bincorr = cmat)  +
                          bindata::rmvbin(n, margprob = meanparam2, bincorr = cmat))
  return(x)
}


