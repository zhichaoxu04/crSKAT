#' genData_pwr.R
#'
#' Generate the data for power simulation studies.
#'
#' @param n the number of observations to generate (sample size).
#' @param alpha1 1st parameter for Gompertz model.
#' @param alpha2 2nd parameter for Gompertz model.
#' @param beta1 3rd parameter for Gompertz model.
#' @param beta2 4th parameter for Gompertz model.
#' @param gMatCausal the genotype matrix of causal SNPs.
#' @param effectSizes the effect size of each causal SNP.
#'
#' @return A list with the elements:
#' \item{tempTime}{n*1 vector of exact event times.}
#' \item{tempType}{n*1 vector of exact event types.}
#'
#' @importFrom stats rbinom
#' @importFrom stats runif
#'
#' @export
#' @examples
#' n <- 1000
#' alpha1 <- -0.058
#' alpha2 <- -0.035
#' beta1 <- 0.03
#' beta2 <- log(1 - exp(beta1 / alpha1)) * alpha2
#' gMatCausal <- gen_GMat(n, q=5, rho=0.1)
#' effectSizes <- seq(0.1, 0.5, length.out=5)
#' genData_pwr(n, alpha1, alpha2, beta1, beta2, gMatCausal, effectSizes, seed=NULL)
#'
#'
genData_pwr <- function(n, alpha1, alpha2, beta1, beta2, gMatCausal, effectSizes, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if(length(effectSizes) != ncol(gMatCausal)){stop("Length of effect Sizes is not equal to the genotype matrix")}
  geneticVec <- gMatCausal %*% effectSizes

  pi1 <- 1 - exp(beta1 / alpha1)
  pi2 <- 1 - exp(beta2 / alpha2)
  tempType <- stats::rbinom(n=n, size=1, prob=pi2) + 1
  tempUnif <- stats::runif(n=n)
  tempTime <- ifelse(tempType == 1, log( 1 - alpha1 / beta1 * log((1 - tempUnif * pi1)) / exp(geneticVec) ) / alpha1,
                     log( 1 - alpha2 / beta2 * log((1 - tempUnif * pi2)) / exp(geneticVec) ) / alpha2)

  return(list(tempTime = tempTime, tempType = tempType))
}
